//! Percent-identity benchmark: profmark-style ROC over Pfam families embedded
//! in a Swissprot decoy background.
//!
//! Construction and analysis live here because they are specific to this
//! benchmark; only execution is shared, via the `run` library.

mod build;
mod parse;

use std::path::PathBuf;

use anyhow::{bail, Result};
use clap::{Parser, Subcommand};

use run::{Asset, Bin, Ctx, Numa, Options, Search};

/// This benchmark's directory, relative to the repository root.


#[derive(Parser)]
#[command(name = "pct-id", about = "percent-identity benchmark")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Assemble a benchmark directory from the profmark split.
    Build(build::Args),
    /// Execute the runs declared in bench.toml.
    Run(RunArgs),
    /// Turn results into the tables the plot scripts consume.
    #[command(subcommand)]
    Parse(parse::Cmd),
}

#[derive(Parser, Debug)]
pub struct RunArgs {
    /// Which benchmark to run, naming benchmark-<size>/.
    #[arg(short, long, default_value = "toy")]
    pub size: String,

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

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Build(args) => build::main(args),
        Command::Run(args) => run_main(args),
        Command::Parse(cmd) => parse::main(cmd),
    }
}

fn run_main(args: RunArgs) -> Result<()> {
    let root = run::repo(env!("CARGO_MANIFEST_DIR"));
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    let config = run::Config::from_path(dir.join("bench.toml"))?;
    let opts = Options {
        filter: args.filter,
        threads: args.threads,
        numa_node: args.numa_node,
        jobs: 1,
        dry_run: args.dry_run,
    };

    let runs = run::plan(&config, &opts)?;
    let bench_dir = dir.join(format!("benchmark-{}", args.size));

    if opts.dry_run {
        // a single search, so listing the runs is the whole matrix
        run::describe(&runs, &[Search::new("", bench_dir.join("target.fa"))]);
        return Ok(());
    }

    if !bench_dir.is_dir() {
        bail!(
            "benchmark directory {} does not exist; run `pct-id build --size {}` first",
            bench_dir.display(),
            args.size
        );
    }

    // absolute, so the cmd column of the runs table can be pasted into a shell
    let bench_dir = bench_dir.canonicalize()?;

    // one target, so no label is needed to disambiguate outputs
    let searches = vec![
        Search::new("", bench_dir.join("target.fa"))
            .with(Asset::Hmm, bench_dir.join("query.hmm"))
            .with(Asset::Sto, bench_dir.join("query.sto"))
            .with(Asset::Fasta, bench_dir.join("query.fa"))
            .with(Asset::Afa, bench_dir.join("afa")),
    ];

    let results = bench_dir.join("results");
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
        bin: Bin::new(root.join("tools/bin")),
        tmp: bench_dir.join("tmp"),
        runs_table: results.join("runs.tbl"),
        results,
        numa,
    };

    run::execute(&config, &runs, &searches, &ctx, opts.jobs)
}
