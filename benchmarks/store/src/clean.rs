//! Removing what a run left behind.
//!
//! Takes the same `paths.toml` a benchmark reads and the same label, and
//! removes what that label names as produced: the run, the analysis and the
//! scratch. The set is left alone -- a build is expensive, and nothing here
//! knows whether it can be made again.

use std::collections::BTreeMap;
use std::path::PathBuf;

use anyhow::bail;
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

/// What a label names, split the one way a clean cares about.
///
/// Every path in a label is something the benchmark produced, except the set,
/// which it reads. That rule needs no vocabulary, which is what lets this
/// clean up after cutoffs -- whose labels name a decoy set and two stage
/// directories rather than the `run` every other benchmark writes.
#[derive(Deserialize, Debug)]
struct Produced {
    set: Option<PathBuf>,

    #[serde(flatten)]
    made: BTreeMap<String, PathBuf>,
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

    // leaked so the names can borrow for the length of the call, which is the
    // whole program: a label's keys are whatever the benchmark chose
    let mut targets: Vec<(&str, PathBuf)> = p
        .made
        .into_iter()
        .map(|(key, path)| (&*key.leak(), file.at(path)))
        .collect();

    match (args.all, p.set) {
        (true, Some(set)) => targets.push(("set", file.at(set))),
        (true, None) => bail!("label {label:?} names no set to remove"),
        (false, _) => {}
    }

    util::clean::run(file.path().parent().unwrap_or(std::path::Path::new(".")), &targets)
}
