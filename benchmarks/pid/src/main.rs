//! Percent-identity benchmark: profmark-style ROC over Pfam families embedded
//! in a Swissprot decoy background.
//!
//! The set is built by `build-set --in pid-toy|pid-real` and read from the
//! store like every other benchmark's. One label fills in every path, so a
//! subcommand never names a directory.

mod inputs;
mod parse;
mod plot;
mod run;
mod search;

use std::path::PathBuf;

use clap::{Parser, Subcommand};
use serde::Deserialize;

#[derive(Parser)]
#[command(name = "pid", about = "percent-identity benchmark")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Search every tool against the benchmark.
    Run(run::Args),
    /// Turn results into the tables the plot scripts consume.
    #[command(subcommand)]
    Parse(parse::Cmd),
    /// Draw the figures from what parse wrote.
    Plot(plot::Args),
}

/// Where this benchmark reads and writes, out of `paths.toml`.
#[derive(Deserialize)]
pub struct Paths {
    /// The set to search.
    pub set: PathBuf,
    /// Where the results and the ledger go.
    pub run: PathBuf,
    /// Where the tables parse works out go.
    pub analysis: PathBuf,
    /// Where the figures go, which is outside the store: a pdf is read by a
    /// person rather than by another pipeline.
    pub figures: PathBuf,
    /// Scratch, and nothing worth keeping.
    pub tmp: PathBuf,
}

impl Paths {
    /// The label named, with every path resolved against the file that named
    /// it.
    pub fn open(label: &str) -> anyhow::Result<Paths> {
        let file = util::paths::File::open(env!("CARGO_MANIFEST_DIR"))?;
        let p: Paths = file.get(label)?;

        Ok(Paths {
            set: file.at(p.set),
            run: file.at(p.run),
            analysis: file.at(p.analysis),
            figures: file.at(p.figures),
            tmp: file.at(p.tmp),
        })
    }

    /// What to print when no label was named.
    pub fn listing(usage: &str) -> anyhow::Result<String> {
        Ok(util::paths::File::open(env!("CARGO_MANIFEST_DIR"))?.listing(usage))
    }
}

const USAGE: &str = "pid <run|parse|plot> --in <label>";

impl Command {
    fn label(&self) -> Option<&str> {
        match self {
            Command::Run(a) => a.label.as_deref(),
            Command::Plot(a) => a.label.as_deref(),
            Command::Parse(parse::Cmd::Recall(a)) => a.label.as_deref(),
            Command::Parse(parse::Cmd::Cells(a)) => a.label.as_deref(),
            Command::Parse(parse::Cmd::Score(a)) => a.label.as_deref(),
            Command::Parse(parse::Cmd::Table(a)) => a.label.as_deref(),
        }
    }
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    let label = match cli.command.label() {
        Some(label) => label.to_string(),
        None => {
            println!("{}", Paths::listing(USAGE)?);
            return Ok(());
        }
    };
    let paths = Paths::open(&label)?;

    match cli.command {
        Command::Run(args) => run::main(args, &paths),
        Command::Parse(cmd) => parse::main(cmd, &paths),
        Command::Plot(args) => plot::main(args, &paths),
    }
}
