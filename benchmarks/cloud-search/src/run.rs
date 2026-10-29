//! What nail's cloud search pruning costs, as -A and -B move.
//!
//! Seeds once against one shard, then searches every (A, B) cell off those
//! seeds. -A and -B do nothing until after seeding, so the seed set is the
//! same in every cell; paying mmseqs once and replaying it into all of them is
//! what keeps the runtime axis about cloud search rather than about the seeder.
//!
//! Every cell runs `--reps` times, with the reps as the outer loop, so load
//! from other users on the box spreads across cells rather than landing on
//! whichever ran that minute. Seeding runs once, ahead of every rep.
//!
//! One shard rather than the whole target set: this benchmark is asking how a
//! knob trades off, and the answer doesn't need more sequences than it takes
//! to see the trade.

use std::path::PathBuf;

use anyhow::{Context, ensure};
use clap::Parser;

use michi::{Cmd, PipelineBuilder, Progress, Step, Table};

use util::ledger;
use util::manifest;
use util::search::{Bins, Dirs, SEED_MODE, SEED_S};
use util::set::Set;

/// What the unpruned run is called, which every cell is measured against.
pub const FULL: &str = "full";

/// The ledger setting holding a run's repetition, from 1.
pub const REP: &str = "rep";

/// What the one seeding is called. Every cell replays it, so the name is the
/// same on every row and only the shard moves.
const SEEDING: &str = "once";

/// The -a every cell runs at: nail's own default, passed rather than left off
/// so the manifest records what ran instead of whatever the binary defaulted
/// to that day.
const ATTEMPTS: u32 = 5;

#[derive(Parser, Debug)]
pub struct Args {
    /// Which label of paths.toml to run under. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// Local score pruning thresholds to sweep (nail's -A)
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "2,3,5,8,10,14,20,30,40",
        value_name = "X,X,..."
    )]
    alpha: Vec<f32>,

    /// Global score pruning thresholds to sweep (nail's -B)
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "4,6,8,12,16,22,32,48,64",
        value_name = "X,X,..."
    )]
    beta: Vec<f32>,

    /// Threads per search, and the cores each search is pinned to
    #[arg(short, long)]
    threads: Option<usize>,

    /// Times each cell is run
    #[arg(long)]
    reps: Option<usize>,

    #[arg(long)]
    tmp: Option<PathBuf>,

    #[arg(long)]
    dry_run: bool,
}

/// One point on the grid: a pair of pruning thresholds, or the unpruned
/// reference.
///
/// `Full` is nail with `--full-dp`, which skips the cloud stage and fills the
/// whole matrix. It is the ceiling every pruned cell is measured against: the
/// most nail can find off a given seed set, and the longest it can take to
/// find it.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Cell {
    Pruned { a: f32, b: f32 },
    Full,
}

impl Cell {
    /// What one rep of this cell's results file is called, and what its step
    /// is named.
    fn label(self, rep: usize) -> String {
        match self {
            Cell::Pruned { a, b } => format!("A{a:.1}-B{b:.1}-r{rep}"),
            Cell::Full => format!("{FULL}-r{rep}"),
        }
    }
}

/// Every A against every B, then the unpruned cell on the end.
///
/// Cells where A >= B are in here and are expected to come out identical to
/// each other: A prunes against the best score on the current anti-diagonal
/// and B against the best score anywhere, so a local threshold above the
/// global one never binds. They are left in as a check rather than skipped as
/// waste.
fn cells(alphas: &[f32], betas: &[f32]) -> Vec<Cell> {
    let mut out: Vec<Cell> = alphas
        .iter()
        .flat_map(|&a| betas.iter().map(move |&b| Cell::Pruned { a, b }))
        .collect();

    out.push(Cell::Full);
    out
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    let threads = args.threads.context("--threads is required")?;
    let reps = args.reps.context("--reps is required")?;
    ensure!(reps > 0, "--reps needs to be at least 1");

    ensure!(!args.alpha.is_empty(), "--alpha needs at least one value");
    ensure!(!args.beta.is_empty(), "--beta needs at least one value");

    let bins = Bins::find()?;

    let mut dirs = Dirs::new(&paths.run, &paths.tmp);
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let set = Set::load_needing(&paths.set, crate::NEEDS)?;
    ensure!(set.units().count() > 0, "{} is empty", paths.set.display());

    let mut pl = PipelineBuilder::new().pool(threads).step(dirs.mkdir());

    let units: Vec<_> = set
        .units()
        .map(|unit| Ok((unit.name().to_string(), unit.query_hmm()?, unit.target()?)))
        .collect::<anyhow::Result<_>>()?;

    // every unit gets the whole grid. On a `fixed` set that is the one shard
    // the label names; on a `cross` it is each target source in turn, which is
    // what makes the surface comparable across kinds of target rather than
    // only across sizes of one
    for (shard, query_hmm, target) in &units {
        pl = pl.step(util::search::seed(
            &bins.nail,
            &bins.mmseqs,
            query_hmm,
            target,
            shard,
            &dirs.seeds(SEEDING, shard),
            &dirs,
            threads,
            &util::search::Seeding::new(SEED_S, SEED_MODE),
            &[
                (manifest::STAGE, util::search::SEED.to_string()),
                (manifest::SEEDS, SEEDING.to_string()),
            ],
        ));
    }

    for rep in 1..=reps {
        for (shard, query_hmm, target) in &units {
            // the cells run in the order the grid gives them, which is ascending,
            // so the cheap corner lands first and a sweep that gets killed still
            // leaves a usable surface behind
            for cell in cells(&args.alpha, &args.beta) {
                let label = cell.label(rep);

                let cmd = Cmd::new(&bins.nail)
                    .sub("search")
                    // nail looks for mmseqs at startup even when it is replaying
                    // seeds and will never call it, and nothing here is on PATH
                    .arg("--mmseqs-path", &bins.mmseqs)
                    .arg("-t", threads)
                    .arg("--seeds", dirs.seeds(SEEDING, shard))
                    .arg("--seeds-format", "names")
                    .arg("-E", util::search::EVALUE)
                    .arg("--tmp-dir", dirs.tmp.join("cell"))
                    .arg("--tbl-out", dirs.table(&label, shard))
                    .flag("--allow-overwrite")
                    .field(manifest::NAME, &label)
                    .field(manifest::TOOL, "nail")
                    // every cell replays the one seeding this unit got
                    .field(manifest::SEEDS, SEEDING)
                    .field(manifest::SHARD, shard)
                    .field(REP, rep);

                let cmd = match cell {
                    // the fields are written the way the label is, so a
                    // whole-numbered threshold keeps its decimal point and
                    // `A=2.0` reads against `A2.0-B4.0` rather than beside it
                    Cell::Pruned { a, b } => cmd
                        .arg("-A", a)
                        .arg("-B", b)
                        .arg("-a", ATTEMPTS)
                        .field("A", format!("{a:.1}"))
                        .field("B", format!("{b:.1}")),
                    Cell::Full => cmd.flag("--full-dp"),
                };

                pl = pl.step(
                    Step::serial([cmd.path(query_hmm).path(target)])
                        .name(&label)
                        // every cell on the same cores, so the only thing moving
                        // between them is -A and -B. without this a cell is timed
                        // against whatever else the scheduler ran that second, and
                        // the differences here are small enough for that to show
                        .cores(threads),
                );
            }
        }
    }

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
