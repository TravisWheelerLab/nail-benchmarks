//! Timing every primitive over the rungs of the ladder.
//!
//! Every command here is built by the same [`crate::search`] helpers the real
//! pipelines use, so what gets timed is what will run rather than a model of
//! it. The one thing this adds is a `part` field on each command, naming which
//! primitive it is, so [`super::fit`] can group the rows without parsing names.
//!
//! The loops run ascending on both axes, so the cheap corner lands first and a
//! sweep that gets killed still leaves a usable surface behind.
//!
//! What is deliberately *not* here: every `(A, B)` cell of cloud-search's grid.
//! A replay's cost is driven by the seed set, and the grid is 81 points on a
//! surface that a handful of samples describe. [`super::fit`] fits that surface
//! and [`super::predict`] sums it over the whole grid.

use std::path::PathBuf;

use anyhow::{Context, ensure};
use clap::Parser;

use michi::{Cmd, PipelineBuilder, Progress, Step, Table};

use crate::inputs;
use crate::search::{self, Bins, Dirs, Split};
use util::ledger;
use util::manifest;

use super::Part;

/// The field naming which primitive a run is.
//
// only runs need it. the four primitives that are stages of
// the pipeline -- split, createdb, seed, convert -- are
// already named by their stage, and Part::key is that name
pub const PART: &str = "part";

/// The field naming the query rung, since the shard is the target rung.
pub const QUERY: &str = "q";

/// The field naming the seeding a replay replays.
//
// said outright rather than parsed back out of the run's
// name: a cell is called A2.0-B4.0, so the name is full of
// dots and there is no separator left to split on
pub const SEEDS: &str = "seeds";

/// The seeding mode the replay pipelines use.
//
// not swept: the mode selects a different seeder rather than
// varying a parameter of one
const SEED_MODE: &str = "prog";

/// What `--max-seqs` recall holds mmseqs to.
const MMSEQS_MAX_SEQS: usize = 2000;

#[derive(Parser, Debug)]
pub struct Args {
    /// Threads per search, and the cores each search is pinned to
    #[arg(short, long, default_value_t = 8)]
    threads: usize,

    /// How many times to time each measurement. The fit takes the median
    #[arg(long, default_value_t = 1)]
    reps: usize,

    /// Seeding sensitivities to measure, as nail's --mmseqs-s. The replay
    /// pipelines seed at 12.0
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "12.0",
        value_name = "X,X,..."
    )]
    seed_s: Vec<String>,

    /// Pruning cells to measure, as A:B pairs, plus the unpruned ceiling.
    /// Samples of a surface, not the grid itself
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "2:4,10:16,40:64",
        value_name = "A:B,..."
    )]
    cells: Vec<String>,

    /// nail's prefilter sensitivities, for recall's standalone column
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "10.0",
        value_name = "X,X,..."
    )]
    nail_s: Vec<String>,

    /// mmseqs' sensitivities, for recall's mmseqs column
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "7.5",
        value_name = "X,X,..."
    )]
    mmseqs_s: Vec<String>,

    /// Where the scratch goes
    #[arg(long)]
    tmp: Option<PathBuf>,

    /// Print the commands without running them
    #[arg(long)]
    dry_run: bool,
}

/// One pruning sample: a pair of thresholds, or the unpruned ceiling.
#[derive(Clone, Copy, Debug)]
enum Cell {
    Pruned { a: f32, b: f32 },
    Full,
}

impl Cell {
    fn parse(text: &str) -> anyhow::Result<Cell> {
        if text == "full" {
            return Ok(Cell::Full);
        }

        let (a, b) = text
            .split_once(':')
            .with_context(|| format!("{text:?} is not an A:B pair"))?;

        Ok(Cell::Pruned {
            a: a.parse()
                .with_context(|| format!("{a:?} is not a threshold"))?,
            b: b.parse()
                .with_context(|| format!("{b:?} is not a threshold"))?,
        })
    }

    fn label(self) -> String {
        match self {
            Cell::Pruned { a, b } => format!("A{a:.1}-B{b:.1}"),
            Cell::Full => "full".to_string(),
        }
    }
}

pub fn main(args: Args) -> anyhow::Result<()> {
    ensure!(
        args.threads.is_multiple_of(search::HMMER_CPU),
        "--threads needs to be a multiple of {} (for hmmer)",
        search::HMMER_CPU
    );
    ensure!(args.reps > 0, "--reps needs to be at least 1");

    let cells: Vec<Cell> = args
        .cells
        .iter()
        .map(|text| Cell::parse(text))
        .collect::<anyhow::Result<_>>()?;
    ensure!(!cells.is_empty(), "--cells needs at least one");

    let bins = Bins::find()?;

    let mut dirs = Dirs::new("calibrate");
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let queries = inputs::ladder::query_rungs()?;
    let targets = inputs::ladder::target_rungs()?;

    // split once per query rung: the parts don't depend on the target, and a
    // rung is searched against every one of them
    let splits: Vec<Split> = queries
        .iter()
        .map(|&q| {
            Split::new(
                inputs::ladder::query_hmm(q),
                dirs.tmp.join(format!("hmmer-query/{q}")),
                search::jobs(args.threads),
            )
        })
        .collect();

    let scratch = dirs.tmp.join("scratch");
    let cell_dir = scratch.join("cell");
    let target_db = scratch.join("targetDB/targetDB");

    let mut pl = PipelineBuilder::new().step(dirs.mkdir());
    for (&q, split) in queries.iter().zip(&splits) {
        // one split per query rung, and none of them searches a shard, so
        // without a stage of its own every rung would land in one ledger row
        pl = pl.step(split.step(&[
            (manifest::STAGE, format!("{}.q{q}", Part::Split)),
            (QUERY, q.to_string()),
        ]));
    }

    for &t in &targets {
        let shard = t.to_string();
        let target_fa = inputs::ladder::target(t);

        pl = pl.step(
            Step::serial([
                Cmd::new("rm").name("clean").flag("-rf").path(&scratch),
                Cmd::new("mkdir")
                    .name("dirs")
                    .flag("-p")
                    .path(scratch.join("targetDB")),
                // read the whole target once so the page-cache
                // cost is not charged to whichever tool runs first
                Cmd::new("cat").name("warm").path(&target_fa),
                search::createdb(&bins.mmseqs, &target_fa, &target_db, &shard, args.threads),
            ])
            .name(format!("prep.t{t}")),
        );

        for (&q, split) in queries.iter().zip(&splits) {
            let query_hmm = inputs::ladder::query_hmm(q);
            let query_db = inputs::ladder::query_db(q);

            for rep in 1..=args.reps {
                // the rep is part of the name, not only a field: two rows of
                // one run against one shard cannot disagree on a param, so a
                // rep that were only a field would not distill
                let at = format!("q{q}.r{rep}");
                let named = |what: &str| format!("{what}.{at}");

                let run = |cmd: Cmd, name: &str, tool: &str, part: Part| {
                    cmd.field(manifest::NAME, name)
                        .field(manifest::TOOL, tool)
                        .field(manifest::SHARD, &shard)
                        .field(PART, part.key())
                        .field(QUERY, q)
                };

                pl = pl.step(
                    Step::serial([
                        // at the front, so a cell that failed leaves its
                        // scratch behind and the next one still starts clean
                        Cmd::new("rm").name("clean").flag("-rf").path(&cell_dir),
                        Cmd::new("mkdir")
                            .name("dirs")
                            .flag("-p")
                            .path(cell_dir.join("nail"))
                            .path(cell_dir.join("mmseqs/alnDB")),
                    ])
                    .name(format!("prep.{at}.t{t}")),
                );

                // ---- seeding, and the replays off each seed set

                for s in &args.seed_s {
                    let label = format!("s{s}.{at}");
                    let seeds = dirs.seeds(&format!("{label}.{t}"));

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
                            s,
                            SEED_MODE,
                            // one seeding per (query rung, sensitivity, rep),
                            // and a stage is keyed by its name and its shard
                            &format!("{}.{label}", search::SEED),
                            &[(QUERY, q.to_string())],
                        )
                        .name(format!("seed.{label}.t{t}")),
                    );

                    for cell in &cells {
                        let name = format!("{}.{label}", cell.label());

                        let cmd = Cmd::new(&bins.nail)
                            .sub("search")
                            // nail looks for mmseqs at startup even when it is
                            // replaying seeds and will never call it, and
                            // nothing here is on PATH
                            .arg("--mmseqs-path", &bins.mmseqs)
                            .arg("-t", args.threads)
                            .arg("--seeds", &seeds)
                            .arg("-E", search::EVALUE)
                            .arg("--tmp-dir", cell_dir.join("replay"))
                            .arg("--tbl-out", dirs.table(&name, &shard))
                            .flag("--allow-overwrite");

                        let cmd = match cell {
                            Cell::Pruned { a, b } => cmd
                                .arg("-A", *a)
                                .arg("-B", *b)
                                .field("A", format!("{a:.1}"))
                                .field("B", format!("{b:.1}")),
                            Cell::Full => cmd.flag("--full-dp"),
                        };

                        pl = pl.step(
                            Step::serial([run(
                                cmd.path(&query_hmm)
                                    .path(&target_fa)
                                    .field("s", s)
                                    .field(SEEDS, &label),
                                &name,
                                "nail",
                                Part::Replay,
                            )])
                            .name(format!("replay.{name}.t{t}"))
                            .cores(args.threads),
                        );
                    }
                }

                // ---- nail end to end, which is recall's column

                for s in &args.nail_s {
                    let name = named(&format!("nail-s{s}"));

                    pl = pl.step(
                        Step::serial([run(
                            Cmd::new(&bins.nail)
                                .sub("search")
                                .arg("--mmseqs-path", &bins.mmseqs)
                                .arg("-t", args.threads)
                                .arg("--tmp-dir", cell_dir.join("nail"))
                                .arg("--mmseqs-s", s)
                                .arg("--seed-mode", SEED_MODE)
                                .arg("-E", search::EVALUE)
                                .arg("--tbl-out", dirs.table(&name, &shard))
                                .flag("--allow-overwrite")
                                .path(&query_hmm)
                                .path(&target_fa)
                                .field("s", s),
                            &name,
                            "nail",
                            Part::Nail,
                        )])
                        .name(format!("nail.s{s}.{at}.t{t}"))
                        .cores(args.threads),
                    );
                }

                // ---- mmseqs, which is recall's other column

                for s in &args.mmseqs_s {
                    let name = named(&format!("mmseqs-s{s}"));
                    let cmds = search::Mmseqs {
                        bin: &bins.mmseqs,
                        query_db: &query_db,
                        target_db: &target_db,
                        aln_db: cell_dir.join("mmseqs/alnDB/alnDB"),
                        work: cell_dir.join("mmseqs/work"),
                        out: dirs.table(&name, &shard),
                        threads: args.threads,
                        s: Some(s.clone()),
                        max_seqs: Some(MMSEQS_MAX_SEQS),
                    }
                    .cmds();

                    pl = pl.step(
                        Step::serial([
                            run(cmds.search.field("s", s), &name, "mmseqs", Part::Mmseqs),
                            // the conversion is a stage: what a column cost is
                            // the search, not the pass that reformats it
                            cmds.convert
                                .field(manifest::STAGE, format!("{}.s{s}.{at}", Part::Convert))
                                .field(manifest::SHARD, &shard)
                                .field(QUERY, q),
                        ])
                        .name(format!("mmseqs.s{s}.{at}.t{t}"))
                        .cores(args.threads),
                    );
                }

                // ---- hmmer

                let name = named("hmmer");
                let hmmer = search::hmmer(
                    &bins.hmmsearch,
                    split,
                    &dirs,
                    &name,
                    &shard,
                    &target_fa,
                    &[
                        (PART, Part::Hmmer.key().to_string()),
                        (QUERY, q.to_string()),
                    ],
                );

                pl = pl
                    .step(hmmer.search.name(format!("hmmer.{at}.t{t}")))
                    .step(hmmer.cat.name(format!("cat.{at}.t{t}")));
            }
        }
    }

    pl = pl.step(Cmd::new("rm").name("clean").flag("-rf").path(&scratch));

    let pipeline = pl
        .stderr_dir(dirs.tmp.join("stderr"))
        .sink(Progress::new())
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
