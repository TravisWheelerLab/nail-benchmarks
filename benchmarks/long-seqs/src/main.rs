//! Long-sequence benchmark: six very large protein pairs, used to measure the
//! fraction of the DP matrix nail computes at extreme sequence lengths.
//!
//! Queries and targets are paired rather than crossed, and the inputs are
//! checked into git, so this benchmark has no build step -- [`inputs`] is where
//! they are, not how they were made.

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

/// Where this benchmark reads and writes, as one label of `paths.toml` names
/// them.
#[derive(serde::Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct Paths {
    pub set: std::path::PathBuf,
    pub run: std::path::PathBuf,
    pub analysis: std::path::PathBuf,
    pub tmp: std::path::PathBuf,
}

impl Paths {
    pub fn open(label: &str) -> anyhow::Result<Paths> {
        let file = util::paths::File::open(env!("CARGO_MANIFEST_DIR"))?;
        let p: Paths = file.get(label)?;

        Ok(Paths {
            set: file.at(p.set),
            run: file.at(p.run),
            analysis: file.at(p.analysis),
            tmp: file.at(p.tmp),
        })
    }

    pub fn listing(usage: &str) -> anyhow::Result<String> {
        Ok(util::paths::File::open(env!("CARGO_MANIFEST_DIR"))?.listing(usage))
    }
}

const USAGE: &str = "long-seqs <run|parse> --in <label>";

fn main() -> anyhow::Result<()> {
    let cmd = Cli::parse().command;

    let label = match &cmd {
        Command::Run(a) => a.label.as_deref(),
        Command::Parse(parse::Cmd::Cells(a)) => a.label.as_deref(),
    };

    let Some(label) = label else {
        println!("{}", Paths::listing(USAGE)?);
        return Ok(());
    };
    let paths = Paths::open(label)?;

    match cmd {
        Command::Run(args) => run::main(args, &paths),
        Command::Parse(cmd) => parse::main(cmd, &paths),
    }
}
