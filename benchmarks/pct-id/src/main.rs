//! Percent-identity benchmark: profmark-style ROC over Pfam families embedded
//! in a Swissprot decoy background.

mod build;
mod parse;
mod run;

use std::path::PathBuf;

use clap::{Parser, Subcommand};

pub fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

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
    /// Search every tool against the benchmark.
    Run(run::Args),
    /// Turn results into the tables the plot scripts consume.
    #[command(subcommand)]
    Parse(parse::Cmd),
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        Command::Build(args) => build::main(args),
        Command::Run(args) => run::main(args),
        Command::Parse(cmd) => parse::main(cmd),
    }
}
