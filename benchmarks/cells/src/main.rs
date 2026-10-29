//! How much of the dynamic-programming matrix nail's cloud search fills,
//! against how long the query and target are, on targets as drawn, written
//! backwards, and shuffled.

mod parse;
mod plot;
mod run;

use std::path::PathBuf;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "cells",
    about = "nail's cell fraction, with and without signal"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Search every form of the targets once with nail.
    Run(run::Args),
    /// Turn each form's table into its cell fractions.
    Parse(parse::Args),
    /// Draw the cell fractions of every form on one figure.
    Plot(plot::Args),
}

/// Where this benchmark reads and writes, as one label of `paths.toml` names
/// them.
#[derive(serde::Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct Paths {
    pub set: PathBuf,
    pub run: PathBuf,
    pub analysis: PathBuf,
    pub figures: PathBuf,
    pub tmp: PathBuf,
}

impl Paths {
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

    pub fn listing(usage: &str) -> anyhow::Result<String> {
        Ok(util::paths::File::open(env!("CARGO_MANIFEST_DIR"))?.listing(usage))
    }
}

const USAGE: &str = "cells <run|parse|plot> --in <label>";

fn main() -> anyhow::Result<()> {
    let cmd = Cli::parse().command;

    let label = match &cmd {
        Command::Run(a) => a.label.as_deref(),
        Command::Parse(a) => a.label.as_deref(),
        Command::Plot(a) => a.label.as_deref(),
    };

    let Some(label) = label else {
        println!("{}", Paths::listing(USAGE)?);
        return Ok(());
    };
    let paths = Paths::open(label)?;

    match cmd {
        Command::Run(args) => run::main(args, &paths),
        Command::Parse(args) => parse::main(args, &paths),
        Command::Plot(args) => plot::main(args, &paths),
    }
}
