//! Removing what a run left behind.
//!
//! Takes the same `paths.toml` a benchmark reads and the same label, and
//! removes what that label names as produced: the run, the analysis and the
//! scratch. The set is left alone -- a build is expensive, and nothing here
//! knows whether it can be made again.

use std::path::PathBuf;

use clap::Parser;
use serde::Deserialize;

#[derive(Parser, Debug)]
pub struct Args {
    /// A benchmark's paths.toml
    #[arg(long, value_name = "paths.toml")]
    paths: PathBuf,

    /// Which label of it to take back. Omit to list them
    #[arg(long = "in", value_name = "label")]
    label: Option<String>,

    /// Also remove the set the label reads, which cost a `build-set` run
    #[arg(long)]
    all: bool,
}

/// The keys this needs of a label. A benchmark's own type may name more; the
/// ones not here are no business of a clean.
#[derive(Deserialize, Debug)]
struct Produced {
    set: PathBuf,
    run: PathBuf,
    analysis: PathBuf,
    tmp: PathBuf,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let dir = args
        .paths
        .parent()
        .unwrap_or(std::path::Path::new("."))
        .to_owned();

    let file = util::paths::File::open(dir)?;

    let Some(label) = args.label else {
        println!("{}", file.listing("store clean --paths <paths.toml> --in <label>"));
        return Ok(());
    };

    let p: Produced = file.get(&label)?;

    let mut targets = vec![
        ("run", file.at(p.run)),
        ("analysis", file.at(p.analysis)),
        ("tmp", file.at(p.tmp)),
    ];

    if args.all {
        targets.push(("set", file.at(p.set)));
    }

    util::clean::run(file.path().parent().unwrap_or(std::path::Path::new(".")), &targets)
}
