mod build;
mod cell;
mod parse;
mod plot;
mod run;
mod scores;

use std::path::PathBuf;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "cloud-search",
    about = "what nail's cloud search pruning costs, as -A and -B move"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Draw the query set and the target set into a benchmark directory.
    Build(build::Args),
    /// Seed once, then search every (A, B) cell off those seeds.
    Run(run::Args),
    /// Turn a finished sweep into tables.
    #[command(subcommand)]
    Parse(parse::Cmd),
    /// Draw the figures from grid.tbl.
    Plot(plot::Args),
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        Command::Build(args) => build::main(args),
        Command::Run(args) => run::main(args),
        Command::Parse(cmd) => parse::main(cmd),
        Command::Plot(args) => plot::main(args),
    }
}

pub fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}
