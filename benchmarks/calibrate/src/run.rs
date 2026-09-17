//! Timing every search over the rungs of the ladder.
//!
//! Every command here is built by the same [`search`] helpers the real
//! pipelines use, so what gets timed is what will run rather than a model of
//! it. The one thing this adds is a `part` field on each command, naming which
//! search it is, so [`super::fit`] can group the rows without parsing names.
//!
//! The query is all of Pfam at every rung; only the target grows. The rungs
//! run ascending, so the cheap end lands first and a sweep that gets killed
//! still leaves a usable line behind.
//!
//! nail is timed in two halves at each sensitivity, seeding and then the
//! alignment off those seeds, because that is what the pipelines are built
//! out of and the two add up to a whole nail search.

use std::collections::BTreeMap;
use std::path::PathBuf;

use anyhow::{Context, ensure};
use clap::Parser;

use michi::{Cmd, PipelineBuilder, Progress, Sink, Step, Table};

use util::set::Set;

use search::sweeps::{MMSEQS_MAX_SEQS, MMSEQS_S, NAIL_S, SEED_MODE};
use search::{self, Bins, Dirs, Split};
use util::ledger;
use util::manifest;

use super::Part;

/// The field naming which search a run is.
pub const PART: &str = "part";

#[derive(Parser, Debug)]
pub struct Args {
    /// Which label of paths.toml to run under. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// Threads per search, and the cores each search is pinned to
    #[arg(short, long, default_value_t = 8)]
    threads: usize,

    /// How many times to time each measurement
    #[arg(long, default_value_t = 1)]
    reps: usize,

    /// nail's --mmseqs-s values, seeded and aligned at each. Defaults to what
    /// recall sweeps
    #[arg(long, value_delimiter = ',', default_value = NAIL_S, value_name = "X,X,...")]
    nail_s: Vec<String>,

    /// mmseqs' -s values. Defaults to what recall sweeps
    #[arg(long, value_delimiter = ',', default_value = MMSEQS_S, value_name = "X,X,...")]
    mmseqs_s: Vec<String>,

    /// Which searches to time, comma separated. Every one by default
    #[arg(long, value_delimiter = ',', value_name = "PART,...")]
    parts: Vec<Part>,

    /// Where the scratch goes
    #[arg(long)]
    tmp: Option<PathBuf>,

    /// Print the commands without running them
    #[arg(long)]
    dry_run: bool,
}

/// What each target rung came to, as the build counted it.
fn rung_residues(units: &[util::set::Unit<'_>]) -> anyhow::Result<BTreeMap<usize, u64>> {
    units
        .iter()
        .map(|u| Ok((u.need("target_rung")?.parse()?, u.number("target_residues")?)))
        .collect()
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    ensure!(
        args.threads.is_multiple_of(search::HMMER_CPU),
        "--threads needs to be a multiple of {} (for hmmer)",
        search::HMMER_CPU
    );
    ensure!(args.reps > 0, "--reps needs to be at least 1");

    let parts = match args.parts.is_empty() {
        true => Part::ALL.to_vec(),
        false => args.parts.clone(),
    };
    let timing = |part: Part| parts.contains(&part);

    // an alignment reads the seed list the seeding beside it wrote, so asking
    // for one without the other leaves it nothing to align
    ensure!(
        !timing(Part::Align) || timing(Part::Seed),
        "--parts align needs seed too: an alignment runs off the seed list \
         that seeding writes"
    );

    let bins = Bins::find()?;

    let mut dirs = Dirs::new(&paths.run, &paths.tmp);
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let set = Set::load_as(&paths.set, crate::SHAPE)?;

    let rungs = set.values("query_rung");
    let [query] = &rungs[..] else {
        anyhow::bail!(
            "the ladder has {} query rungs; calibrate searches with all of Pfam at \
             every target size, so it wants exactly one: rebuild with \
             `mgy build ladder --queries 20795`",
            rungs.len()
        );
    };

    // one query rung, so its units are the target axis and nothing else moves
    let units: Vec<_> = set.where_attr("query_rung", query).collect();
    anyhow::ensure!(!units.is_empty(), "no units at query rung {query}");

    let query_hmm = units[0].query_hmm()?;
    let query_db = units[0].query_db()?;

    // one split for the whole sweep: the parts do not depend on the target,
    // and hmmer searches every rung with the same ones
    let split = Split::new(
        &query_hmm,
        dirs.tmp.join("hmmer-query"),
        search::jobs(args.threads),
    );

    let scratch = dirs.tmp.join("scratch");
    let target_db = scratch.join("targetDB/targetDB");

    let mut pl = PipelineBuilder::new().step(dirs.mkdir());
    if timing(Part::Hmmer) {
        pl = pl.step(split.step(&[]));
    }

    for unit in &units {
        let shard = unit.need("target_rung")?.to_string();
        let target_fa = unit.target()?;

        let mut prep = vec![
            Cmd::new("rm").name("clean").flag("-rf").path(&scratch),
            // one alnDB per sensitivity: mmseqs refuses to write over a
            // database an earlier search left, so sharing one would fail the
            // second run of every rung
            args.mmseqs_s.iter().fold(
                Cmd::new("mkdir")
                    .name("dirs")
                    .flag("-p")
                    .path(scratch.join("targetDB")),
                |cmd, s| cmd.path(scratch.join(format!("mmseqs/alnDB-s{s}"))),
            ),
            // read the whole target once so the page-cache cost is not charged
            // to whichever tool runs first
            Cmd::new("cat").name("warm").path(&target_fa),
        ];

        if timing(Part::Mmseqs) {
            prep.push(search::createdb(
                &bins.mmseqs,
                &target_fa,
                &target_db,
                &shard,
                args.threads,
            ));
        }

        pl = pl.step(Step::serial(prep).name(format!("prep.t{shard}")));

        for rep in 1..=args.reps {
            let run = |cmd: Cmd, name: &str, tool: &str, part: Part| {
                cmd.field(manifest::NAME, name)
                    .field(manifest::TOOL, tool)
                    .field(manifest::SHARD, &shard)
                    .field(PART, part.key())
            };

            // ---- nail, three timings per sensitivity

            for s in &args.nail_s {
                let at = format!("s{s}.r{rep}");
                let seeds = dirs.seeds(&format!("{at}.{shard}"));

                let name = format!("seed-{at}");
                if timing(Part::Seed) {
                    pl = pl.step(
                        search::seed(
                            &bins.nail,
                            &bins.mmseqs,
                            &query_hmm,
                            &target_fa,
                            &shard,
                            &seeds,
                            &dirs,
                            args.threads,
                            &search::Seeding::new(s, SEED_MODE),
                            &[
                                (manifest::NAME, name.clone()),
                                (manifest::TOOL, "nail".to_string()),
                                (PART, Part::Seed.key().to_string()),
                                ("s", s.clone()),
                            ],
                        )
                        .name(format!("{name}.t{shard}"))
                        .cores(args.threads),
                    );
                }

                let name = format!("align-{at}");
                if timing(Part::Align) {
                    pl = pl.step(
                        Step::serial([run(
                            Cmd::new(&bins.nail)
                                .sub("search")
                                // nail looks for mmseqs at startup even when it
                                // is aligning a seed set and will never call
                                // it, and nothing here is on PATH
                                .arg("--mmseqs-path", &bins.mmseqs)
                                .arg("-t", args.threads)
                                .arg("--seeds", &seeds)
                                .arg("-E", search::EVALUE)
                                .arg("--tmp-dir", scratch.join("align"))
                                .arg("--tbl-out", dirs.table(&name, &shard))
                                .flag("--allow-overwrite")
                                .path(&query_hmm)
                                .path(&target_fa)
                                .field("s", s),
                            &name,
                            "nail",
                            Part::Align,
                        )])
                        .name(format!("{name}.t{shard}"))
                        .cores(args.threads),
                    );
                }
            }

            // ---- mmseqs

            for s in args.mmseqs_s.iter().filter(|_| timing(Part::Mmseqs)) {
                let name = format!("mmseqs-s{s}.r{rep}");
                let cmds = search::Mmseqs {
                    bin: &bins.mmseqs,
                    query_db: &query_db,
                    target_db: &target_db,
                    aln_db: scratch.join(format!("mmseqs/alnDB-s{s}/alnDB")),
                    work: scratch.join(format!("mmseqs/work-s{s}")),
                    out: dirs.table(&name, &shard),
                    threads: args.threads,
                    s: Some(s.clone()),
                    max_seqs: Some(MMSEQS_MAX_SEQS),
                }
                .cmds();

                pl = pl.step(
                    Step::serial([
                        run(cmds.search.field("s", s), &name, "mmseqs", Part::Mmseqs),
                        // the conversion carries no part, so it is timed into
                        // the manifest and left out of the model: what a
                        // column cost is the search
                        cmds.convert.field(manifest::SHARD, &shard),
                    ])
                    .name(format!("{name}.t{shard}"))
                    .cores(args.threads),
                );
            }

            // ---- hmmer

            if !timing(Part::Hmmer) {
                continue;
            }

            let name = format!("hmmer.r{rep}");
            let hmmer = search::hmmer(
                &bins.hmmsearch,
                &split,
                &dirs,
                &name,
                &shard,
                &target_fa,
                &[(PART, Part::Hmmer.key().to_string())],
            );

            pl = pl
                .step(hmmer.search.name(format!("{name}.t{shard}")))
                .step(hmmer.cat.name(format!("cat.{name}.t{shard}")));
        }
    }

    pl = pl.step(Cmd::new("rm").name("clean").flag("-rf").path(&scratch));

    let pipeline = pl
        .stderr_dir(dirs.tmp.join("stderr"))
        .sink(Progress::new())
        .sink(Rungs::new(rung_residues(&units)?))
        .sink(Table::new(dirs.root.join("manifest.tbl")))
        .build()
        .context("failed to build the sweep")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    // the ledger describes the results this run is about to replace, so it
    // goes before the run rather than after the failure of one
    ledger::clear(&dirs.root);
    pipeline.run()?;
    ledger::record(&dirs.root)
}

/// Prints each rung's timings once that rung is done, rather than leaving a
/// sweep of several hours silent until the ledger lands.
struct Rungs {
    residues: BTreeMap<usize, u64>,
    at: Option<String>,
    rows: Vec<(String, f64)>,
}

impl Rungs {
    fn new(residues: BTreeMap<usize, u64>) -> Rungs {
        Rungs {
            residues,
            at: None,
            rows: Vec::new(),
        }
    }

    fn flush(&mut self) {
        let Some(rung) = self.at.take() else { return };
        let rows = std::mem::take(&mut self.rows);

        let residues = rung
            .parse()
            .ok()
            .and_then(|rung| self.residues.get(&rung).copied())
            .unwrap_or(0);

        let total: f64 = rows.iter().map(|(_, wall)| wall).sum();
        let width = rows.iter().map(|(name, _)| name.len()).max().unwrap_or(0);

        println!("\n  rung {rung} -- {residues} target residues");
        for (name, wall) in &rows {
            println!("    {name:<width$}  {wall:>9.2} s");
        }
        println!("    {:<width$}  {total:>9.2} s\n", "rung total");
    }
}

impl Sink for Rungs {
    fn step_done(&mut self, step: &Step) -> anyhow::Result<()> {
        // a step with no part is setup -- the scratch, the warm read, the
        // database, the cat -- and is not one of the timings being reported
        let Some(item) = step.items().next() else {
            return Ok(());
        };

        let fields = item.fields();
        let (Some(rung), Some(name)) = (fields.get(manifest::SHARD), fields.get(manifest::NAME))
        else {
            return Ok(());
        };
        if !fields.contains_key(PART) {
            return Ok(());
        }

        if self.at.as_deref() != Some(rung.as_str()) {
            self.flush();
            self.at = Some(rung.clone());
        }

        // the step's own wall clock rather than the command's, so hmmer's
        // parts count once between them the way the ledger folds them
        if let Some(wall) = step.wall_s() {
            self.rows.push((name.clone(), wall));
        }

        Ok(())
    }

    fn finish(&mut self) -> anyhow::Result<()> {
        self.flush();
        Ok(())
    }
}
