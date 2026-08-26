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
//! -a is left alone. It is a real knob and it carries most of the weight at
//! aggressive thresholds, but it is not what this benchmark is asking about.

use std::path::PathBuf;

use anyhow::{Context, ensure};
use clap::Parser;

use pipeline::{Cmd, PipelineBuilder, Progress, Step, Table};

use crate::manifest;
use crate::search::{self, Bins, Dirs, Split};

/// The seeding settings this benchmark searches with.
const MMSEQS_S: &str = "12.0";
const SEED_MODE: &str = "prog";

#[derive(Parser, Debug)]
pub struct Args {
    /// Which target shard to search.
    #[arg(long, default_value = "1", value_name = "N")]
    shard: String,

    /// Local score pruning thresholds to sweep (nail's -A).
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "2,3,5,8,10,14,20,30,40",
        value_name = "X,X,..."
    )]
    alpha: Vec<f32>,

    /// Global score pruning thresholds to sweep (nail's -B).
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "4,6,8,12,16,22,32,48,64",
        value_name = "X,X,..."
    )]
    beta: Vec<f32>,

    /// Threads per search, and the cores each search is pinned to.
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
    Pruned { a: f32, b: f32 },
    Full,
}

impl Cell {
    /// What this cell's results file is called, and what its step is named.
    fn label(self) -> String {
        match self {
            Cell::Pruned { a, b } => format!("A{a:.1}-B{b:.1}"),
            Cell::Full => "full".to_string(),
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

pub fn main(args: Args) -> anyhow::Result<()> {
    ensure!(
        args.threads.is_multiple_of(search::HMMER_CPU),
        "--threads needs to be a multiple of {} (for hmmer)",
        search::HMMER_CPU
    );
    ensure!(!args.alpha.is_empty(), "--alpha needs at least one value");
    ensure!(!args.beta.is_empty(), "--beta needs at least one value");

    let bins = Bins::find()?;

    let mut dirs = Dirs::new("cloud-search");
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let query_hmm = crate::queries().join("query.hmm");
    let target = search::shard(&args.shard);

    ensure!(
        query_hmm.is_file() && target.is_file(),
        "{} or {} is missing; run `mgy build` first",
        query_hmm.display(),
        target.display()
    );

    let split = Split::new(
        &query_hmm,
        dirs.tmp.join("hmmer-query"),
        search::jobs(args.threads),
    );

    let mut pl = PipelineBuilder::new()
        .step(dirs.mkdir())
        .step(split.step())
        .step(search::seed(
            &bins.nail,
            &bins.mmseqs,
            &query_hmm,
            &target,
            &args.shard,
            &dirs,
            args.threads,
            MMSEQS_S,
            SEED_MODE,
        ));

    for step in search::truth(&bins.hmmsearch, &split, &dirs, &args.shard, &target) {
        pl = pl.step(step);
    }

    // the cells run in the order the grid gives them, which is ascending, so
    // the cheap corner lands first and a sweep that gets killed still leaves a
    // usable surface behind
    for cell in cells(&args.alpha, &args.beta) {
        let label = cell.label();

        let cmd = Cmd::new(&bins.nail)
            .sub("search")
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
            Cell::Pruned { a, b } => cmd
                .arg("-A", a)
                .arg("-B", b)
                .field("A", format!("{a:.1}"))
                .field("B", format!("{b:.1}")),
            Cell::Full => cmd.flag("--full-dp"),
        };

        pl = pl.step(
            Step::serial([cmd.path(&query_hmm).path(&target)])
                .name(&label)
                // every cell on the same cores, so the only thing moving
                // between them is -A and -B. without this a cell is timed
                // against whatever the scheduler felt like that second, and
                // the differences here are small enough for that to show
                .cores(args.threads),
        );
    }

    let pipeline = pl
        .stderr_dir(dirs.tmp.join("stderr"))
        .sink(Progress::new())
        .sink(Table::new(dirs.root.join("runs.tbl")))
        .build()
        .context("failed to build the sweep")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    pipeline.run()
}
