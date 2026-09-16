//! How much of what hmmer finds nail and mmseqs find, as their prefilter
//! sensitivity moves.
//!
//! Reads a [`shape::FIXED`] set: one query set against target shards of equal
//! size. A shard is a unit of work rather than a variable, so the sweep spans
//! every one of them and the query never moves.

mod run;

use std::path::PathBuf;

use clap::{Parser, Subcommand};

/// The set this reads. Named here rather than inside the run, so what this
/// benchmark needs of a dataset is the first thing in the file.
pub const SHAPE: &util::set::Shape = &util::set::shape::FIXED;

#[derive(Parser)]
#[command(name = "recall", about = "recall against prefilter sensitivity")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Search nail and mmseqs against every shard, sweeping their prefilter
    /// sensitivity.
    Run(run::Args),
    /// Turn the results into a table, and that into numbers.
    #[command(subcommand)]
    Parse(Parse),
}

/// The analyses this benchmark's results answer.
///
/// `scores.tbl` is one row per pair with a column per tool, which is the honest
/// shape for a prefilter sweep: moving a prefilter changes which pairs a tool
/// reports, not what it scores them.
#[derive(Subcommand)]
enum Parse {
    /// Read the results into scores.tbl, one row per pair, one score column
    /// per tool.
    Scores(scores::parse::ScoresArgs),
    /// What every run found and what it cost, one row per run.
    Summary(scores::parse::TableArgs),
}

/// Where this benchmark reads and writes, as one label of `paths.toml` names
/// them.
#[derive(serde::Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct Paths {
    /// The set to search.
    pub set: PathBuf,
    /// Where the results and the ledger go.
    pub run: PathBuf,
    /// Where the tables parse works out go.
    pub analysis: PathBuf,
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
            tmp: file.at(p.tmp),
        })
    }

    /// What to print when no label was named.
    pub fn listing(usage: &str) -> anyhow::Result<String> {
        Ok(util::paths::File::open(env!("CARGO_MANIFEST_DIR"))?.listing(usage))
    }
}

const USAGE: &str = "recall <run|parse> --in <label>";

impl Command {
    fn label(&self) -> Option<&str> {
        match self {
            Command::Run(a) => a.label.as_deref(),
            Command::Parse(Parse::Scores(a)) => a.label.as_deref(),
            Command::Parse(Parse::Summary(a)) => a.label.as_deref(),
        }
    }

    fn run(self, paths: &Paths) -> anyhow::Result<()> {
        match self {
            Command::Run(args) => run::main(args, paths),
            Command::Parse(Parse::Scores(mut a)) => {
                a.run = paths.run.clone();
                a.set = paths.set.clone();
                a.analysis = paths.analysis.clone();
                scores::parse::main(scores::parse::Cmd::Scores(a))
            }
            Command::Parse(Parse::Summary(mut a)) => {
                a.analysis = paths.analysis.clone();
                scores::parse::main(scores::parse::Cmd::Summary(a))
            }
        }
    }
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    // one label fills in every path, so a subcommand never names a directory
    let label = match cli.command.label() {
        Some(label) => label.to_string(),
        None => {
            println!("{}", Paths::listing(USAGE)?);
            return Ok(());
        }
    };
    let paths = Paths::open(&label)?;

    cli.command.run(&paths)
}
