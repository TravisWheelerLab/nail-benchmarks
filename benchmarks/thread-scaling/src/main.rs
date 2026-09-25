//! Thread scaling: how nail, mmseqs and hmmer speed up as they are given more
//! cores, over a fixed amount of work.
//!
//! Strong scaling, so the set never moves: every unit of a [`shape::CROSS`]
//! set is searched at every rung of a thread ladder, by each tool, a few times
//! over. What comes out is each tool's wall clock against its thread count,
//! and whether the hits it reports hold still as the count changes.
//!
//! [`shape::CROSS`]: util::set::shape::CROSS

mod parse;
mod plot;
mod run;

use std::path::PathBuf;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "thread-scaling", about = "wall clock against thread count")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Search every unit at every thread count, with every tool.
    Run(run::Args),
    /// Write scaling.tbl and agree.tbl from a finished run.
    Parse(parse::Args),
    /// Draw speedup and efficiency against threads from scaling.tbl.
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

/// This crate's directory, fixed at compile time. The plot script hangs off it.
pub fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

const USAGE: &str = "thread-scaling <run|parse|plot> --in <label>";

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
