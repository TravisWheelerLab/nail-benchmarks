mod build;
mod cutoffs;
mod parse;
mod run;
mod util;

use anyhow::{bail, Context};
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "mgnify", about = "mgnify/pfam benchmark")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Shard the MGnify proteins into a benchmark directory.
    Build(build::Args),
    /// Execute the runs declared in bench.toml against a range of shards.
    Run(run::Args),
    /// Learn per-family false-positive score cutoffs from reversed decoys.
    #[command(subcommand)]
    Cutoffs(cutoffs::Cmd),
    /// Turn results into tables and analyses.
    #[command(subcommand)]
    Parse(parse::Cmd),
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        Command::Build(args) => build::main(args),
        Command::Run(args) => run::main(args),
        Command::Cutoffs(cmd) => cutoffs::main(cmd),
        Command::Parse(cmd) => parse::main(cmd),
    }
}
