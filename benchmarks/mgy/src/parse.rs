//! Turning a finished pipeline into a table, and that table into the numbers.
//!
//! `scores` and `summary` are recall's. The shard-parallel collector under
//! them is not -- any pipeline's results go through it -- but what a row looks
//! like is a decision each analysis makes for itself, so this is where the
//! pipelines stop being interchangeable.

use std::path::{Path, PathBuf};

use anyhow::{Context, bail};
use clap::{Parser, Subcommand};

use crate::analyze;
use crate::inputs;
use crate::scores;

#[derive(Subcommand)]
pub enum Cmd {
    /// Read recall's results into scores.tbl, one row per pair, one score
    /// column per tool.
    Scores(ScoresArgs),
    /// Read a sweep's results into runs.tbl, one row per pair, one score
    /// column per run, and whether seeding found the pair.
    Runs(ScoresArgs),
    /// What every run found and what it cost, one row per run.
    Summary(TableArgs),
    /// Where the hits hmmer found were lost, one row per run per
    /// checkpoint.
    Funnel(TableArgs),
}

pub fn main(cmd: Cmd) -> anyhow::Result<()> {
    match cmd {
        Cmd::Scores(args) => scores(args),
        Cmd::Runs(args) => runs(args),
        Cmd::Summary(args) => summary(args),
        Cmd::Funnel(args) => funnel(args),
    }
}

#[derive(Parser, Debug)]
pub struct ScoresArgs {
    /// A pipeline directory, or the name of one under benchmarks/mgy/outputs/
    #[arg(value_name = "pipeline")]
    pipeline: String,

    /// The per-family cutoffs a calibration learned. Defaults to the committed
    /// ones
    #[arg(long, value_name = "cutoffs.tbl")]
    cutoffs: Option<PathBuf>,

    /// Which decoy to cut at. The cutoffs file holds each family's five
    /// best-scoring decoys, so this admits at most `c` false positives per
    /// family. It is fixed here rather than in the analyses, since the cutoff
    /// travels in the table
    #[arg(short = 'c', default_value_t = 2, value_name = "N")]
    c: usize,

    /// The query set that was searched. Defaults to the shared one
    #[arg(long, value_name = "query.hmm")]
    queries: Option<PathBuf>,

    /// The directory holding the target shards. Defaults to the shared one
    #[arg(long, value_name = "dir")]
    targets: Option<PathBuf>,

    #[arg(short, long, value_name = "scores.tbl")]
    out: Option<PathBuf>,

    /// How many shards to collect at once. Defaults to the machine's cores
    #[arg(long, value_name = "N")]
    threads: Option<usize>,

    /// How many gigabytes the collectors may hold between them. Defaults to
    /// half of what the machine has
    #[arg(long, value_name = "GB")]
    mem: Option<f64>,
}

fn scores(args: ScoresArgs) -> anyhow::Result<()> {
    let set = Inputs::resolve(args, "scores.tbl")?;

    let count = scores::write::collect(scores::write::Args {
        dir: &set.dir,
        query_hmm: &set.query_hmm,
        targets: &set.targets,
        cutoffs: &set.cutoffs,
        c: set.c,
        out: &set.out,
        threads: set.threads,
        mem: set.mem,
    })?;

    set.report(count);
    Ok(())
}

fn runs(args: ScoresArgs) -> anyhow::Result<()> {
    let set = Inputs::resolve(args, "runs.tbl")?;

    let count = scores::runs::collect(scores::runs::Args {
        dir: &set.dir,
        query_hmm: &set.query_hmm,
        targets: &set.targets,
        cutoffs: &set.cutoffs,
        c: set.c,
        out: &set.out,
        threads: set.threads,
        mem: set.mem,
    })?;

    set.report(count);
    Ok(())
}

/// What both collectors are pointed at, once the defaults are filled in.
struct Inputs {
    dir: PathBuf,
    query_hmm: PathBuf,
    targets: PathBuf,
    cutoffs: PathBuf,
    c: usize,
    out: PathBuf,
    threads: usize,
    mem: u64,
}

impl Inputs {
    fn resolve(args: ScoresArgs, name: &str) -> anyhow::Result<Inputs> {
        let dir = pipeline(&args.pipeline)?;

        let cutoffs = match args.cutoffs {
            Some(path) => path,
            None => util::tools::mgy_cutoffs()?,
        };

        let threads = match args.threads {
            Some(threads) => threads,
            None => std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1),
        };

        let mem = match args.mem {
            Some(gb) => (gb * (1u64 << 30) as f64) as u64,
            None => scores::collect::ram() / 2,
        };

        Ok(Inputs {
            out: args.out.unwrap_or_else(|| dir.join(name)),
            dir,
            query_hmm: args.queries.unwrap_or_else(inputs::fixed::query_hmm),
            targets: args.targets.unwrap_or_else(inputs::fixed::targets),
            cutoffs,
            c: args.c,
            threads,
            mem,
        })
    }

    fn report(&self, count: scores::shard::Count) {
        println!(
            "wrote {} ({} rows out of {} hits)",
            self.out.display(),
            count.rows,
            count.hits
        );
    }
}

#[derive(Parser, Debug)]
pub struct TableArgs {
    /// The scores.tbl `parse scores` wrote, the pipeline directory holding
    /// one, or the name of one under benchmarks/mgy/outputs/
    #[arg(value_name = "recall|cloud-search|hit-loss")]
    scores: String,

    #[arg(short, long, value_name = "out.tbl")]
    out: Option<PathBuf>,
}

fn summary(args: TableArgs) -> anyhow::Result<()> {
    // either grammar: a summary is the one analysis both pipelines' tables
    // answer, so naming a directory has to find whichever it wrote
    let path = table(&args, &["scores.tbl", "runs.tbl"])?;
    let out = beside(&path, args.out, "summary.tbl")?;

    analyze::summary(&path, &out)?;

    println!("wrote {}", out.display());
    Ok(())
}

fn funnel(args: TableArgs) -> anyhow::Result<()> {
    let path = table(&args, &["runs.tbl"])?;
    let out = beside(&path, args.out, "funnel.tbl")?;

    analyze::funnel(&path, &out)?;

    println!("wrote {}", out.display());
    Ok(())
}

/// Where the table an analysis reads is, given either as a path or by
/// pipeline.
///
/// `names` are the tables this analysis can read, in the order it prefers
/// them. A path is taken as given; a pipeline is searched for the first of
/// them it holds.
fn table(args: &TableArgs, names: &[&str]) -> anyhow::Result<PathBuf> {
    let given = PathBuf::from(&args.scores);
    if given.is_file() {
        return Ok(given);
    }

    let dir = pipeline(&args.scores)?;

    match names.iter().map(|name| dir.join(name)).find(|p| p.is_file()) {
        Some(path) => Ok(path),
        None => {
            let wanted: Vec<String> = names
                .iter()
                .map(|name| {
                    let (stem, _) = name.split_once('.').unwrap_or((name, ""));
                    format!("`mgy parse {stem}`")
                })
                .collect();

            bail!(
                "no {} in {}; run {} first",
                names.join(" or "),
                dir.display(),
                wanted.join(" or "),
            )
        }
    }
}

fn beside(path: &Path, out: Option<PathBuf>, name: &str) -> anyhow::Result<PathBuf> {
    match out {
        Some(path) => Ok(path),
        None => Ok(path
            .parent()
            .context("the table has no directory")?
            .join(name)),
    }
}

/// A pipeline directory, given either as a path or by name.
pub fn pipeline(name: &str) -> anyhow::Result<PathBuf> {
    let dir = match PathBuf::from(name) {
        path if path.is_dir() => path,
        _ => crate::outputs().join(name),
    };

    if !dir.is_dir() {
        bail!("no pipeline directory at {}", dir.display());
    }

    Ok(dir)
}
