//! MGnify/Pfam benchmark: the full Pfam profile set searched against
//! metagenomic proteins, sharded so the work can be split across NUMA nodes.
//!
//! Construction and analysis live here; only execution is shared, via the
//! `run` library.

mod build;
mod fasta;
mod parse;

use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};

use run::{Asset, Bin, Ctx, Numa, Options, Search};

use build::DIR;

#[derive(Parser)]
#[command(name = "mgnify", about = "mgnify/pfam benchmark")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Shard the MGnify proteins, optionally writing reversed decoy shards.
    Build(build::Args),
    /// Execute the runs declared in bench.toml against a range of shards.
    Run(RunArgs),
    /// Turn results into tables and learned cutoffs.
    #[command(subcommand)]
    Parse(parse::Cmd),
}

#[derive(Parser, Debug)]
pub struct RunArgs {
    /// Which benchmark to run, naming benchmark-<size>/.
    #[arg(short, long, default_value = "toy")]
    pub size: String,

    /// Shard range, as `N` or `FIRST-LAST` (1-based, inclusive). Defaults to
    /// every shard in the benchmark. The two-node production runs split this,
    /// e.g. 1-500 on one node and 501-1000 on the other.
    #[arg(long)]
    pub shards: Option<String>,

    /// Search the reversed decoy shards instead of the real ones.
    #[arg(long)]
    pub rev: bool,

    /// Only run entries whose name matches this glob.
    #[arg(short, long)]
    pub filter: Option<String>,

    /// Override the thread count from bench.toml.
    #[arg(short, long)]
    pub threads: Option<usize>,

    /// Pin to a NUMA node. Absent means no pinning and no numactl call.
    #[arg(long)]
    pub numa_node: Option<usize>,

    /// Where results are written.
    #[arg(long)]
    pub results: Option<PathBuf>,

    /// Scratch directory.
    #[arg(long)]
    pub tmp: Option<PathBuf>,

    /// List the expanded runs and exit without executing anything.
    #[arg(long)]
    pub dry_run: bool,

    /// Repository root, if it cannot be discovered automatically.
    #[arg(long)]
    pub root: Option<PathBuf>,
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Build(args) => build::main(args),
        Command::Run(args) => run_main(args),
        Command::Parse(cmd) => parse::main(cmd),
    }
}

fn count_shards(dir: &std::path::Path) -> usize {
    std::fs::read_dir(dir)
        .map(|entries| {
            entries
                .filter_map(|e| e.ok())
                .filter(|e| e.path().extension().is_some_and(|x| x == "fa"))
                .count()
        })
        .unwrap_or(0)
}

/// Parse `N` or `FIRST-LAST` into an inclusive range.
fn shard_range(spec: &str) -> Result<(usize, usize)> {
    let (a, b) = match spec.split_once('-') {
        Some((a, b)) => (a.trim(), b.trim()),
        None => (spec.trim(), spec.trim()),
    };

    let a: usize = a
        .parse()
        .with_context(|| format!("bad shard range {spec:?}"))?;
    let b: usize = b
        .parse()
        .with_context(|| format!("bad shard range {spec:?}"))?;

    if a == 0 {
        bail!("shards are 1-based; {spec:?} starts at 0");
    }
    if a > b {
        bail!("shard range {spec:?} is empty");
    }

    Ok((a, b))
}

fn run_main(args: RunArgs) -> Result<()> {
    let repo = run::repo_root(args.root.as_deref())?;
    let dir = repo.join(DIR);

    let config = run::Config::from_path(dir.join("bench.toml"))?;
    let opts = Options {
        filter: args.filter,
        threads: args.threads,
        numa_node: args.numa_node,
        dry_run: args.dry_run,
    };
    let runs = run::plan(&config, &opts)?;

    let bench = dir.join(format!("benchmark-{}", args.size));
    if !bench.is_dir() {
        bail!(
            "benchmark directory {} does not exist; run `mgnify build --size {}` first",
            bench.display(),
            args.size
        );
    }
    let bench = bench.canonicalize()?;

    let shard_dir = bench.join(if args.rev { "mgy-rev" } else { "mgy" });
    let available = count_shards(&shard_dir);
    if available == 0 {
        bail!(
            "no shards in {}{}",
            shard_dir.display(),
            if args.rev {
                "; rebuild with --reverse"
            } else {
                ""
            }
        );
    }

    let (first, last) = match &args.shards {
        Some(spec) => shard_range(spec)?,
        None => (1, available),
    };
    if last > available {
        bail!("asked for shard {last} but the benchmark has only {available}");
    }

    // the query set is the same for every shard; mmseqs converts it to a
    // profile db once and reuses it, keyed on the source path
    let searches: Vec<Search> = (first..=last)
        .map(|i| {
            Search::new(i.to_string(), shard_dir.join(format!("{i}.fa")))
                .with(Asset::Hmm, bench.join("query.hmm"))
                .with(Asset::Sto, bench.join("query.sto"))
        })
        .collect();

    if opts.dry_run {
        run::describe(&runs, &searches);
        return Ok(());
    }

    let results = args.results.unwrap_or_else(|| bench.join("results"));
    let tmp = args.tmp.unwrap_or_else(|| bench.join("tmp"));
    if results.exists() {
        std::fs::remove_dir_all(&results)?;
    }
    std::fs::create_dir_all(&results)?;

    let threads = runs.iter().map(|r| r.threads).max().unwrap_or(1);
    let numa = match args.numa_node.or(config.defaults.numa_node) {
        Some(node) => Some(Numa::new(node, threads)?),
        None => None,
    };

    let ctx = Ctx {
        bin: Bin::new(repo.join("tools/bin")),
        tmp,
        results,
        numa,
    };

    run::execute(&config, &runs, &searches, &ctx)
}
