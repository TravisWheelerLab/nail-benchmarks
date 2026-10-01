//! Turning a finished pipeline into a table, and that table into the numbers.
//!
//! `scores` and `summary` are recall's. The shard-parallel collector under
//! them is not -- any pipeline's results go through it -- but what a row looks
//! like is a decision each analysis makes for itself, so this is where the
//! pipelines stop being interchangeable.

use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};
use clap::{Parser, Subcommand};

use util::set::Set;

use crate::scores::analyze;

#[derive(Subcommand)]
pub enum Cmd {
    /// Read a sweep's results into runs.tbl, one row per pair, one score
    /// column per run, and whether seeding found the pair.
    Runs(ScoresArgs),
    /// Where the ceiling's hits were lost, one row per run per checkpoint.
    Stages(TableArgs),
    /// How far down the prefilter list the hits sit, one row per run per
    /// depth bin.
    Depth(DepthArgs),
}

pub fn main(cmd: Cmd) -> anyhow::Result<()> {
    match cmd {
        Cmd::Runs(args) => runs(args),
        Cmd::Stages(args) => stages(args),
        Cmd::Depth(args) => depth(args),
    }
}

#[derive(Parser, Debug)]
pub struct ScoresArgs {
    /// Which label of paths.toml to read. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// The run directory to read: `ledger.tbl` and `results/`. Filled in from
    /// the label by the benchmark that owns it.
    #[clap(skip)]
    pub run: PathBuf,

    /// The set that was searched.
    #[clap(skip)]
    pub set: PathBuf,

    /// Where the table goes.
    #[clap(skip)]
    pub analysis: PathBuf,

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

fn runs(args: ScoresArgs) -> anyhow::Result<()> {
    let set = Inputs::resolve(args, "runs.tbl")?;

    let count = crate::scores::runs::collect(crate::scores::runs::Args {
        dir: &set.dir,
        query_hmm: &set.query_hmm,
        set: &set.set,
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
    set: Set,
    query_hmm: PathBuf,
    cutoffs: PathBuf,
    c: usize,
    out: PathBuf,
    threads: usize,
    mem: u64,
}

impl Inputs {
    fn resolve(args: ScoresArgs, name: &str) -> anyhow::Result<Inputs> {
        let dir = args.run;

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
            None => crate::scores::collect::ram() / 2,
        };

        // held to the one column this opens rather than to a shape: the
        // pipelines that write these tables search a `fixed` shard and a
        // `cross` unit alike, and a name in the manifest is not what decides
        // whether the table can be read
        let set = Set::load_needing(&args.set, &[util::set::Rep::QueryHmm])?;

        // the query is a property of the set rather than of a unit, and every
        // pipeline that writes one of these tables searches one query set. A
        // cross may pair several, so that is checked rather than assumed --
        // taking the first would report one query's scores under another's
        let queries: std::collections::BTreeSet<_> = set
            .units()
            .map(|u| u.query_hmm())
            .collect::<anyhow::Result<_>>()?;

        ensure!(
            queries.len() == 1,
            "{} names {} query sets; this reads one",
            args.set.display(),
            queries.len()
        );

        let query_hmm = queries
            .into_iter()
            .next()
            .context("the set names no units")?;

        let out = match args.out {
            Some(path) => path,
            None => args.analysis.join(name),
        };

        // the collectors stream into this rather than building it in memory,
        // so the directory has to exist before the first write rather than
        // after the last
        if let Some(parent) = out.parent() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("failed to create {}", parent.display()))?;
        }

        Ok(Inputs {
            out,
            dir,
            set,
            query_hmm,
            cutoffs,
            c: args.c,
            threads,
            mem,
        })
    }

    fn report(&self, count: crate::scores::shard::Count) {
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
    /// Which label of paths.toml to read. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// The analysis directory holding the table, filled in from the label.
    #[clap(skip)]
    pub analysis: PathBuf,

    #[arg(short, long, value_name = "out.tbl")]
    out: Option<PathBuf>,
}

#[derive(Parser, Debug)]
pub struct DepthArgs {
    /// Which label of paths.toml to read. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// The run directory holding results/, filled in from the label.
    #[clap(skip)]
    pub run: PathBuf,

    /// The analysis directory holding runs.tbl, filled in from the label.
    #[clap(skip)]
    pub analysis: PathBuf,

    /// Ranks in the first bin; each bin after it doubles the ceiling
    #[arg(long, default_value_t = 200, value_name = "N")]
    bin: u64,

    #[arg(short, long, value_name = "out.tbl")]
    out: Option<PathBuf>,
}

fn depth(args: DepthArgs) -> anyhow::Result<()> {
    let table = TableArgs {
        label: args.label,
        analysis: args.analysis,
        out: None,
    };
    let path = table_at(&table, &["runs.tbl"])?;
    let out = beside(&path, args.out, "depth.tbl")?;
    let hits = beside(&out, None, "hits.tbl")?;
    let lists = beside(&out, None, "lists.tbl")?;

    crate::scores::depth::depth(&path, &args.run, args.bin, &out, &hits, &lists)?;

    println!("wrote {}", out.display());
    println!("wrote {}", hits.display());
    println!("wrote {}", lists.display());
    Ok(())
}

fn stages(args: TableArgs) -> anyhow::Result<()> {
    let path = table_at(&args, &["runs.tbl"])?;
    let out = beside(&path, args.out, "stages.tbl")?;

    analyze::stages(&path, &out)?;

    println!("wrote {}", out.display());
    Ok(())
}

/// Where the table an analysis reads is, given either as a path or by
/// pipeline.
///
/// `names` are the tables this analysis can read, in the order it prefers
/// them. A path is taken as given; a pipeline is searched for the first of
/// them it holds.
fn table_at(args: &TableArgs, names: &[&str]) -> anyhow::Result<PathBuf> {
    if args.analysis.is_file() {
        return Ok(args.analysis.clone());
    }

    let dir = &args.analysis;

    match names
        .iter()
        .map(|name| dir.join(name))
        .find(|p| p.is_file())
    {
        Some(path) => Ok(path),
        None => {
            let wanted: Vec<String> = names
                .iter()
                .map(|name| {
                    let (stem, _) = name.split_once('.').unwrap_or((name, ""));
                    format!("`loss-decomp parse {stem}`")
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
