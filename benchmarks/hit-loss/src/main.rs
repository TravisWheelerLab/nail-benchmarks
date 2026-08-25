mod build;
mod parse;
mod run;

use std::path::PathBuf;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "hit-loss",
    about = "which stage of the nail pipeline drops the hits hmmer finds"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Draw the query set and the target set into a benchmark directory.
    Build(build::Args),
    /// Seed once, run hmmer for truth, then run nail once at its defaults.
    Run(run::Args),
    /// Turn a finished run into a funnel table.
    Parse(parse::Args),
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        Command::Build(args) => build::main(args),
        Command::Run(args) => run::main(args),
        Command::Parse(args) => parse::main(args),
    }
}

pub fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}
