//! Every unit searched at every rung of a thread ladder, by every arm, a few
//! times over.
//!
//! Four arms. nail at `-t N`, mmseqs at `--threads N`, and hmmer twice: one
//! `hmmsearch --cpu N`, and the query split N/2 ways with each part at
//! `--cpu 2`, which is how every other benchmark here runs it. Each search is
//! pinned to N cores, so a rung is the thread count and the cores together.
//!
//! The repetitions are the outer loop. Load from other users on the box comes
//! and goes over minutes, and a rep-major order spreads whatever there was
//! across every rung rather than piling it onto whichever rung was running.
//! `/proc/loadavg` is read either side of every search into `load.tbl`, so a
//! point that was run under load can be seen to have been.

use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use anyhow::{Context, ensure};
use clap::Parser;

use michi::{Closure, Cmd, PipelineBuilder, Progress, Step, Table};

use util::ledger;
use util::manifest;
use util::search::{Bins, Dirs, EVALUE, HMMER_CPU, SEED_MODE, SEED_S, Split};
use util::set::Set;
use util::split::Kind;

/// What `load.tbl` is called, beside the ledger.
pub const LOAD: &str = "load.tbl";

/// The ledger setting that names the arm, since two of them are hmmer.
pub const ARM: &str = "arm";

/// The ledger setting holding a run's thread count.
pub const THREADS: &str = "threads";

/// The ledger setting holding a run's repetition, from 1.
pub const REP: &str = "rep";

/// One way of running a tool at N threads.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Arm {
    Nail,
    Mmseqs,
    Hmmer,
    HmmerSplit,
}

impl Arm {
    pub fn name(self) -> &'static str {
        match self {
            Arm::Nail => "nail",
            Arm::Mmseqs => "mmseqs",
            Arm::Hmmer => "hmmer",
            Arm::HmmerSplit => "hmmer-split",
        }
    }

    fn tool(self) -> &'static str {
        match self {
            Arm::Nail => "nail",
            Arm::Mmseqs => "mmseqs",
            Arm::Hmmer | Arm::HmmerSplit => "hmmer",
        }
    }

    /// Whether this arm has a point at `threads`.
    fn runs_at(self, threads: usize) -> bool {
        match self {
            // N/2 parts at --cpu 2 is no parts at all on one thread
            Arm::HmmerSplit => threads >= HMMER_CPU,
            _ => true,
        }
    }
}

/// What a run's table and ledger row are filed under.
pub fn run_name(arm: Arm, threads: usize, rep: usize) -> String {
    format!("{}-t{threads}-r{rep}", arm.name())
}

#[derive(Parser, Debug)]
pub struct Args {
    /// Which label of paths.toml to run under. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// Thread counts to run each arm at, and the cores each search is pinned to
    #[arg(long, value_delimiter = ',', value_name = "N,N,...")]
    rungs: Vec<usize>,

    /// Times each point is run
    #[arg(long, default_value_t = 1)]
    reps: usize,

    /// Arms to run
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "nail,mmseqs,hmmer,hmmer-split",
        value_name = "ARM,ARM,..."
    )]
    arms: Vec<Arm>,

    #[arg(long)]
    tmp: Option<PathBuf>,

    #[arg(long)]
    dry_run: bool,
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    ensure!(!args.rungs.is_empty(), "--rungs needs at least one value");
    ensure!(
        !args.rungs.contains(&0),
        "--rungs holds a 0; a search needs a thread"
    );
    ensure!(args.reps > 0, "--reps needs to be at least 1");

    let mut rungs = args.rungs.clone();
    rungs.sort_unstable();
    rungs.dedup();
    let top = *rungs.last().expect("checked non-empty above");

    let bins = Bins::find()?;

    let mut dirs = Dirs::new(&paths.run, &paths.tmp);
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let set = Set::load_as(&paths.set, &util::set::shape::CROSS)?;
    let units: Vec<_> = set.units().collect();
    ensure!(!units.is_empty(), "{} is empty", paths.set.display());

    let load = Load::new(dirs.root.join(LOAD));

    let mut pl = PipelineBuilder::new().pool(top).step(dirs.clean());

    // ---- what every rep reuses: mmseqs' target databases and hmmer's splits

    let mut prep = Vec::new();
    for unit in &units {
        let shard = unit.name().to_string();
        let scratch = dirs.tmp.join(&shard);
        let target_db = scratch.join("targetDB/targetDB");

        pl = pl.step(
            Step::serial([
                Cmd::new("mkdir")
                    .name("dirs")
                    .flag("-p")
                    .path(scratch.join("targetDB")),
                util::search::createdb(&bins.mmseqs, &unit.target()?, &target_db, &shard, top),
            ])
            .name(format!("prep.{shard}")),
        );

        // one split per rung, since the part count is what the rung moves.
        // the parts do not depend on the rep, so they are cut once
        let mut splits = Vec::new();
        for &threads in rungs.iter().filter(|&&n| Arm::HmmerSplit.runs_at(n)) {
            if !args.arms.contains(&Arm::HmmerSplit) {
                break;
            }
            let split = Split::new(
                unit.query_hmm()?,
                Kind::Hmm,
                scratch.join(format!("split-t{threads}")),
                util::search::jobs(threads),
            );
            pl = pl.step(split.step(
                &format!("split.{shard}.t{threads}"),
                &[(manifest::SHARD, shard.clone())],
            ));
            splits.push((threads, split));
        }

        prep.push((target_db, splits));
    }

    // ---- the ladder

    for rep in 1..=args.reps {
        for (unit, (target_db, splits)) in units.iter().zip(&prep) {
            let shard = unit.name().to_string();
            let query_hmm = unit.query_hmm()?;
            let target = unit.target()?;
            let scratch = dirs.tmp.join(&shard);

            for &threads in &rungs {
                for &arm in args.arms.iter().filter(|arm| arm.runs_at(threads)) {
                    let name = run_name(arm, threads, rep);
                    let fields = [
                        (ARM, arm.name().to_string()),
                        (THREADS, threads.to_string()),
                        (REP, rep.to_string()),
                        (manifest::SHARD, shard.clone()),
                    ];

                    pl = pl.step(load.step(&name, &shard, "start"));

                    match arm {
                        Arm::Nail => {
                            pl = pl.step(
                                util::search::nail(
                                    &bins.nail,
                                    &bins.mmseqs,
                                    &query_hmm,
                                    &target,
                                    &dirs.table(&name, &shard),
                                    &scratch.join(&name),
                                    threads,
                                    EVALUE,
                                    &[
                                        ("--mmseqs-s", SEED_S.to_string()),
                                        ("--seed-mode", SEED_MODE.to_string()),
                                    ],
                                    &[],
                                    &name,
                                    &fields,
                                )
                                .name(format!("{name}.{shard}"))
                                .cores(threads),
                            );
                        }
                        Arm::Mmseqs => {
                            let work = scratch.join(&name);
                            let cmds = util::search::Mmseqs {
                                bin: &bins.mmseqs,
                                query_db: &unit.query_db()?,
                                target_db,
                                aln_db: work.join("alnDB/alnDB"),
                                work: work.join("work"),
                                out: dirs.table(&name, &shard),
                                threads,
                                s: Some(SEED_S.to_string()),
                                max_seqs: None,
                                evalue: EVALUE,
                            }
                            .cmds();

                            pl = pl.step(
                                Step::serial([
                                    Cmd::new("mkdir")
                                        .name("dirs")
                                        .flag("-p")
                                        .path(work.join("alnDB")),
                                    util::search::tag(cmds.search, &name, arm.tool(), &fields),
                                    cmds.convert.field(manifest::SHARD, &shard),
                                    // an alignment db per run, 3 reps x 8
                                    // rungs x 2 units of them, is disk
                                    // nothing reads once the table is out
                                    Cmd::new("rm").name("rm").flag("-rf").path(&work),
                                ])
                                .name(format!("{name}.{shard}"))
                                .cores(threads),
                            );
                        }
                        Arm::Hmmer => {
                            let cmd = Cmd::new(&bins.hmmsearch)
                                .name("hmmsearch")
                                .arg("--cpu", threads)
                                .arg("--tblout", dirs.table(&name, &shard))
                                .arg("-E", EVALUE)
                                .path(&query_hmm)
                                .path(&target);

                            pl = pl.step(
                                Step::serial([util::search::tag(cmd, &name, arm.tool(), &fields)])
                                    .name(format!("{name}.{shard}"))
                                    .cores(threads),
                            );
                        }
                        Arm::HmmerSplit => {
                            let split = &splits
                                .iter()
                                .find(|(n, _)| *n == threads)
                                .expect("a split is cut for every rung this arm runs at")
                                .1;

                            // the builder tags shard itself
                            let hmmer = util::search::hmmer(
                                &bins.hmmsearch,
                                split,
                                &dirs,
                                &name,
                                arm.tool(),
                                &shard,
                                &target,
                                EVALUE,
                                false,
                                &fields[..3],
                            );
                            pl = pl.step(hmmer.search).step(hmmer.cat);
                        }
                    }

                    pl = pl.step(load.step(&name, &shard, "end"));
                }
            }
        }
    }

    let pipeline = pl
        .stderr_dir(dirs.tmp.join("stderr"))
        .sink(Progress::new())
        .sink(Table::new(dirs.root.join("manifest.tbl")))
        .build()
        .context("failed to build the ladder")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    // the ledger and the load log describe the results this run is about to
    // replace, so they go before the run rather than after the failure of one
    ledger::clear(&dirs.root);
    std::fs::remove_file(dirs.root.join(LOAD)).ok();
    pipeline.run()?;
    ledger::record(&dirs.root)
}

/// `/proc/loadavg` at either end of every search, rewritten whole after each
/// reading so a run that dies partway leaves what it saw.
#[derive(Clone)]
struct Load {
    path: PathBuf,
    rows: Arc<Mutex<Vec<[String; 6]>>>,
}

impl Load {
    fn new(path: PathBuf) -> Load {
        Load {
            path,
            rows: Arc::new(Mutex::new(Vec::new())),
        }
    }

    /// A step that reads the load and records it against one run.
    ///
    /// It carries no `name` and no `stage`, so the ledger has no row for it.
    fn step(&self, name: &str, shard: &str, when: &str) -> Step<'static> {
        let load = self.clone();
        let (name, shard, when) = (name.to_string(), shard.to_string(), when.to_string());

        Step::from_closures([Closure::new("load", move || {
            let text =
                std::fs::read_to_string("/proc/loadavg").context("failed to read /proc/loadavg")?;
            let mut avg = text.split_whitespace().map(str::to_string);
            let mut next = || avg.next().unwrap_or_default();

            let mut rows = load.rows.lock().expect("no reading panics holding it");
            rows.push([name, shard, when, next(), next(), next()]);
            Ok(write_load(&load.path, &rows)?)
        })])
        .name("load")
    }
}

fn write_load(path: &Path, rows: &[[String; 6]]) -> anyhow::Result<()> {
    let mut table = toil::Table::new(toil::Schema::new([
        "name", "shard", "when", "load1", "load5", "load15",
    ]));
    for row in rows {
        table.row(row.clone());
    }

    table
        .write(path)
        .with_context(|| format!("failed to write {}", path.display()))
}
