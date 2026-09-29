//! Where nail's wall clock goes, stage by stage, and how that changes with
//! the input. One nail search per unit with `-s` on, and its printed stage
//! tree kept beside the hit table.

mod parse;
mod run;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "stage-times",
    about = "nail's stage times, relative to each other"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Search every unit once, keeping nail's stage tree.
    Run(run::Args),
    /// Turn the stage trees into stages.tbl and counts.tbl.
    Parse(parse::Args),
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

const USAGE: &str = "stage-times <run|parse> --in <label>";

fn main() -> anyhow::Result<()> {
    let cmd = Cli::parse().command;

    let label = match &cmd {
        Command::Run(a) => a.label.as_deref(),
        Command::Parse(a) => a.label.as_deref(),
    };

    let Some(label) = label else {
        println!("{}", Paths::listing(USAGE)?);
        return Ok(());
    };
    let paths = Paths::open(label)?;

    match cmd {
        Command::Run(args) => run::main(args, &paths),
        Command::Parse(args) => parse::main(args, &paths),
    }
}
