//! What nail's cloud search pruning costs, as -A and -B move.
//!
//! Seeds once against one shard, then searches every (A, B) cell off those
//! seeds. -A and -B do nothing until after seeding, so the seed set is the
//! same in every cell; paying mmseqs once and replaying it into all of them is
//! what keeps the runtime axis about cloud search rather than about the seeder.
//!
//! One shard rather than the whole target set: this benchmark is asking how a
//! knob trades off, and the answer doesn't need more sequences than it takes
//! to see the trade.
//!
//! Every grid is run twice, at -a 5 and at -a 0. -a decides how many times
//! nail may retry a pair whose clouds came apart, so a cell that prunes hard
//! can buy back sensitivity somewhere -A and -B do not show. The second arm is
//! what separates the two, and it is not optional: a surface at one -a alone
//! cannot say which of them it is measuring.

use std::path::PathBuf;

use anyhow::{Context, ensure};
use clap::Parser;

use michi::{Cmd, PipelineBuilder, Progress, Step, Table};

use search::sweeps::{SEED_MODE, SEED_S};
use search::{self, Bins, Dirs, Split};
use util::ledger;
use util::manifest;
use util::set::Set;

/// The column hmmer's run becomes, which every cell is measured against.
const HMMER: &str = "hmmer";

/// The -a the whole grid is run at, in order.
///
/// 5 is nail's own default, passed rather than left off so the manifest
/// records what ran instead of whatever the binary defaulted to that day. 0
/// turns the recovery off.
const ATTEMPTS: [u32; 2] = [5, 0];

#[derive(Parser, Debug)]
pub struct Args {
    /// Which label of paths.toml to run under. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// Which unit of the set to search
    #[arg(long, default_value = "1", value_name = "N")]
    shard: String,

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
    #[arg(short, long, default_value_t = 8)]
    threads: usize,

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
    Pruned { a: f32, b: f32, attempts: u32 },
    Full,
}

impl Cell {
    /// What this cell's results file is called, and what its step is named.
    fn label(self) -> String {
        match self {
            Cell::Pruned { a, b, attempts } => format!("A{a:.1}-B{b:.1}-a{attempts}"),
            Cell::Full => "full".to_string(),
        }
    }
}

/// Every A against every B at each -a, then the unpruned cell on the end.
///
/// -a is the outer loop, so each of its values gets a whole grid before the
/// next one starts and a sweep that is killed partway leaves one complete
/// surface rather than a piece of each. The default arm goes first, so what a
/// killed sweep leaves is the surface that is already understood.
///
/// The unpruned cell is emitted once. `--full-dp` fills the matrix and never
/// runs the cloud stage, so there is no disjoint cloud for -a to recover and
/// no second timing to take.
///
/// Cells where A >= B are in here and are expected to come out identical to
/// each other: A prunes against the best score on the current anti-diagonal
/// and B against the best score anywhere, so a local threshold above the
/// global one never binds. They are left in as a check rather than skipped as
/// waste.
fn cells(alphas: &[f32], betas: &[f32]) -> Vec<Cell> {
    let mut out: Vec<Cell> = ATTEMPTS
        .iter()
        .flat_map(|&attempts| {
            alphas.iter().flat_map(move |&a| {
                betas.iter().map(move |&b| Cell::Pruned { a, b, attempts })
            })
        })
        .collect();

    out.push(Cell::Full);
    out
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    ensure!(
        args.threads.is_multiple_of(search::HMMER_CPU),
        "--threads needs to be a multiple of {} (for hmmer)",
        search::HMMER_CPU
    );
    ensure!(!args.alpha.is_empty(), "--alpha needs at least one value");
    ensure!(!args.beta.is_empty(), "--beta needs at least one value");

    let bins = Bins::find()?;

    let mut dirs = Dirs::new(&paths.run, &paths.tmp);
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let set = Set::load_as(&paths.set, crate::SHAPE)?;

    let unit = set
        .units()
        .find(|u| u.name() == args.shard)
        .with_context(|| {
            format!(
                "{} has no unit {:?}; it holds {}",
                paths.set.display(),
                args.shard,
                set.units().map(|u| u.name()).collect::<Vec<_>>().join(", ")
            )
        })?;

    let query_hmm = unit.query_hmm()?;
    let target = unit.target()?;

    let split = Split::new(
        &query_hmm,
        dirs.tmp.join("hmmer-query"),
        search::jobs(args.threads),
    );

    let mut pl = PipelineBuilder::new()
        .step(dirs.mkdir())
        .step(split.step(&[]))
        .step(search::seed(
            &bins.nail,
            &bins.mmseqs,
            &query_hmm,
            &target,
            &args.shard,
            &dirs.seeds(&args.shard),
            &dirs,
            args.threads,
            &search::Seeding::new(SEED_S, SEED_MODE),
            &[(manifest::STAGE, search::SEED.to_string())],
        ));

    let hmmer = search::hmmer(
        &bins.hmmsearch,
        &split,
        &dirs,
        HMMER,
        &args.shard,
        &target,
        &[],
    );
    pl = pl.step(hmmer.search).step(hmmer.cat);

    // the cells run in the order the grid gives them, which is ascending, so
    // the cheap corner lands first and a sweep that gets killed still leaves a
    // usable surface behind
    for cell in cells(&args.alpha, &args.beta) {
        let label = cell.label();

        let cmd = Cmd::new(&bins.nail)
            .sub("search")
            // nail looks for mmseqs at startup even when it is replaying seeds
            // and will never call it, and nothing here is on PATH
            .arg("--mmseqs-path", &bins.mmseqs)
            .arg("-t", args.threads)
            .arg("--seeds", dirs.seeds(&args.shard))
            .arg("-E", search::EVALUE)
            .arg("--tmp-dir", dirs.tmp.join("cell"))
            .arg("--tbl-out", dirs.table(&label, &args.shard))
            .flag("--allow-overwrite")
            .field(manifest::NAME, &label)
            .field(manifest::TOOL, "nail")
            .field(manifest::SHARD, &args.shard);

        let cmd = match cell {
            // the fields are written the way the label is, so a whole-numbered
            // threshold keeps its decimal point and `A=2.0` reads against
            // `A2.0-B4.0` rather than beside it
            Cell::Pruned { a, b, attempts } => cmd
                .arg("-A", a)
                .arg("-B", b)
                .arg("-a", attempts)
                .field("A", format!("{a:.1}"))
                .field("B", format!("{b:.1}"))
                // spelled out rather than `a`, which a reader of the table
                // would have to tell from `A` by its case alone
                .field("attempts", attempts),
            Cell::Full => cmd.flag("--full-dp"),
        };

        pl = pl.step(
            Step::serial([cmd.path(&query_hmm).path(&target)])
                .name(&label)
                // every cell on the same cores, so the only thing moving
                // between them is -A and -B. without this a cell is timed
                // against whatever else the scheduler ran that second, and
                // the differences here are small enough for that to show
                .cores(args.threads),
        );
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
