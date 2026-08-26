//! Long-sequence benchmark: six very large protein pairs, used to measure the
//! fraction of the DP matrix nail computes at extreme sequence lengths.
//!
//! Queries and targets are paired rather than crossed, and the inputs are
//! checked into git, so this benchmark has no build step.

mod parse;

use std::path::PathBuf;

use anyhow::bail;
use clap::{Parser, Subcommand};
use pipeline::{Cmd, PipelineBuilder, Progress, Step, Table};
use tools::nail;

/// This benchmark's directory, fixed at compile time.
pub fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// How many query/target pairs are checked in.
pub const PAIRS: usize = 6;

/// The name every pair's hit table and runs.tbl row is filed under.
pub const RUN_NAME: &str = "nail";

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
    #[arg(short, long, default_value_t = 24)]
    pub threads: usize,

    #[arg(long)]
    pub results: Option<PathBuf>,

    #[arg(long)]
    pub tmp: Option<PathBuf>,

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
    let dir = dir();

    let results = args.results.unwrap_or_else(|| dir.join("results"));
    let tmp = args.tmp.unwrap_or_else(|| dir.join("tmp"));

    let nail_bin = nail()?;

    let mut pl = PipelineBuilder::new().step(Cmd::new("mkdir").flag("-p").path(&results));

    for i in 1..=PAIRS {
        let query = dir.join(format!("query/{i}.query.fa"));
        let target = dir.join(format!("target/{i}.target.fa"));

        if !target.exists() {
            bail!("missing input {}", target.display());
        }

        pl = pl.step(
            Step::serial([Cmd::new(&nail_bin)
                .sub("search")
                .arg("-t", args.threads)
                .arg("--tmp-dir", tmp.join(i.to_string()))
                .flag("--allow-overwrite")
                // widens the sparse band enough that these pairs align at all
                .arg("--f32-p", 5)
                .arg("--tbl-out", results.join(format!("{RUN_NAME}.{i}.tbl")))
                .path(&query)
                .path(&target)
                .field("name", RUN_NAME)
                .field("search", i)])
            .name(i.to_string()),
        );
    }

    let pipeline = pl
        .stderr_dir(tmp.join("stderr"))
        .sink(Progress::new())
        .sink(Table::new(results.join("runs.tbl")))
        .build()?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    pipeline.run()
}
