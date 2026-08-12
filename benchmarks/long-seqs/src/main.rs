//! Long-sequence benchmark: six very large protein pairs, used to measure the
//! fraction of the DP matrix nail computes at extreme sequence lengths.
//!
//! Queries and targets are paired rather than crossed, and the inputs are
//! checked into git, so this benchmark has no build step.

mod parse;

use std::path::PathBuf;

use anyhow::bail;
use clap::{Parser, Subcommand};

use run::{Asset, Bin, Ctx, Numa, Options, Search};

/// This benchmark's directory, fixed at compile time.
pub fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// How many query/target pairs are checked in.
pub const PAIRS: usize = 6;

#[derive(Parser)]
#[command(name = "long-seqs", about = "long sequence benchmark")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Search each query against its paired target.
    Run(RunArgs),
    /// Turn results into the tables the plot scripts consume.
    #[command(subcommand)]
    Parse(parse::Cmd),
}

#[derive(Parser, Debug)]
pub struct RunArgs {
    /// Only run entries whose name matches this glob.
    #[arg(short, long)]
    pub filter: Option<String>,

    /// Override the thread count from bench.toml.
    #[arg(short, long)]
    pub threads: Option<usize>,

    /// Pin to a NUMA node. Absent means no pinning and no numactl call.
    #[arg(long)]
    pub numa_node: Option<usize>,

    /// List the expanded runs and exit without executing anything.
    #[arg(long)]
    pub dry_run: bool,

}

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        Command::Run(args) => run_main(args),
        Command::Parse(cmd) => parse::main(cmd),
    }
}

fn run_main(args: RunArgs) -> anyhow::Result<()> {
    let repo = run::repo(env!("CARGO_MANIFEST_DIR"));
    let dir = dir();

    let config = run::Config::from_path(dir.join("bench.toml"))?;
    let opts = Options {
        filter: args.filter,
        threads: args.threads,
        numa_node: args.numa_node,
        jobs: 1,
        dry_run: args.dry_run,
    };
    let runs = run::plan(&config, &opts)?;

    // zipped, not crossed: pair i is query i against target i
    let searches: Vec<Search> = (1..=PAIRS)
        .map(|i| {
            Search::new(
                i.to_string(),
                dir.join(format!("target/{i}.target.fa")),
            )
            .with(Asset::Fasta, dir.join(format!("query/{i}.query.fa")))
        })
        .collect();

    if opts.dry_run {
        run::describe(&runs, &searches);
        return Ok(());
    }

    for search in &searches {
        if !search.target.exists() {
            bail!("missing input {}", search.target.display());
        }
    }

    let results = dir.join("results");
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
        tmp: dir.join("tmp"),
        runs_table: results.join("runs.tbl"),
        results,
        numa,
    };

    run::execute(&config, &runs, &searches, &ctx, opts.jobs)
}
