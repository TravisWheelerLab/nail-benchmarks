//! What can be done to the store itself, rather than to any one benchmark.
//!
//! Both of these are about runs rather than sets, and neither names a
//! benchmark: `import` writes a ledger for results produced somewhere else,
//! and `clean` takes the store's generated trees back.

mod clean;
mod import;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "store", about = "what a run left behind, and taking it back")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Bring in result tables produced somewhere other than here, timed by
    /// the `.time` files that came back with them.
    Import(import::Args),
    /// Remove the runs, the analyses and the scratch.
    Clean(clean::Args),
}

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        Command::Import(args) => import::main(args),
        Command::Clean(args) => clean::main(args),
    }
}
