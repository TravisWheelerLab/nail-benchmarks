//! Long-sequence benchmark: six very large protein pairs, used to measure the
//! fraction of the DP matrix nail computes at extreme sequence lengths.
//!
//! Queries and targets are paired rather than crossed, and the inputs are
//! checked into git, so this benchmark has no build step -- [`inputs`] is where
//! they are, not how they were made.

mod inputs;
mod parse;
mod run;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "long-seqs", about = "long sequence benchmark")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Search each query against its paired target.
    Run(run::Args),
    /// Turn results into the tables the plot scripts consume.
    #[command(subcommand)]
    Parse(parse::Cmd),
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        Command::Run(args) => run::main(args),
        Command::Parse(cmd) => parse::main(cmd),
    }
}
