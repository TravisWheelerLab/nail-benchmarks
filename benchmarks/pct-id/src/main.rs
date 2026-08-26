//! Percent-identity benchmark: profmark-style ROC over Pfam families embedded
//! in a Swissprot decoy background.

mod build;
mod inputs;
mod parse;
mod plot;
mod run;
mod search;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "pct-id", about = "percent-identity benchmark")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Assemble an input set from the profmark split.
    Build(build::Args),
    /// Search every tool against the benchmark.
    Run(run::Args),
    /// Turn results into the tables the plot scripts consume.
    #[command(subcommand)]
    Parse(parse::Cmd),
    /// Draw the figures from what parse wrote.
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
