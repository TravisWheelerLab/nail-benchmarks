//! What nail's cloud-search pruning costs and what it loses.
//!
//! Reads a [`shape::FIXED`] set and searches one unit of it, seeding once and
//! then searching every `(A, B)` pruning cell off those same seeds, so the
//! pruning parameters are the only thing moving.

mod plot;
mod run;

use std::path::PathBuf;

use clap::{Parser, Subcommand};

/// The set this reads.
pub const SHAPE: &util::set::Shape = &util::set::shape::FIXED;

#[derive(Parser)]
#[command(name = "cloud-search", about = "the (A, B) pruning surface")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Seed once, then search every (A, B) cell off those seeds.
    Run(run::Args),
    /// Turn the results into a table, and that into numbers.
    #[command(subcommand)]
    Parse(Parse),
    /// Draw the pruning heatmaps and the tradeoff curve from summary.tbl.
    Plot(plot::Args),
}

/// `-A` and `-B` constrain the dynamic programming, so two cells can score one
/// pair differently. That is what makes a column per run the honest shape here,
/// where recall wants a column per tool.
#[derive(Subcommand)]
enum Parse {
    /// Read the results into runs.tbl, one row per pair, one score column per
    /// run, and whether seeding found the pair.
    Runs(scores::parse::ScoresArgs),
    /// What every run found and what it cost, one row per run.
    Summary(scores::parse::TableArgs),
    /// Where the hits hmmer found were lost, one row per run per checkpoint.
    Funnel(scores::parse::TableArgs),
}

/// This crate's directory, fixed at compile time. The plot script hangs off
/// it; everything a run produces lives in the store.
pub fn dir() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
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

const USAGE: &str = "cloud-search <run|parse|plot> --in <label>";

impl Command {
    fn label(&self) -> Option<&str> {
        match self {
            Command::Run(a) => a.label.as_deref(),
            Command::Plot(a) => a.label.as_deref(),
            Command::Parse(Parse::Runs(a)) => a.label.as_deref(),
            Command::Parse(Parse::Summary(a)) => a.label.as_deref(),
            Command::Parse(Parse::Funnel(a)) => a.label.as_deref(),
        }
    }

    fn run(self, paths: &Paths) -> anyhow::Result<()> {
        match self {
            Command::Run(args) => run::main(args, paths),
            Command::Plot(args) => plot::main(args, paths),
            Command::Parse(Parse::Runs(mut a)) => {
                a.run = paths.run.clone();
                a.set = paths.set.clone();
                a.analysis = paths.analysis.clone();
                scores::parse::main(scores::parse::Cmd::Runs(a))
            }
            Command::Parse(Parse::Summary(mut a)) => {
                a.analysis = paths.analysis.clone();
                scores::parse::main(scores::parse::Cmd::Summary(a))
            }
            Command::Parse(Parse::Funnel(mut a)) => {
                a.analysis = paths.analysis.clone();
                scores::parse::main(scores::parse::Cmd::Funnel(a))
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
