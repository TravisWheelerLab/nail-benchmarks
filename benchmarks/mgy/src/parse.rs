//! Turning any finished pipeline into the one table.
//!
//! There is nothing per-pipeline in here. What the columns are comes out of
//! `runs.tbl`, so `parse recall` and `parse cloud-search` are the same code
//! pointed at different directories.

use std::path::PathBuf;

use anyhow::bail;
use clap::Parser;

use crate::scores::{Inputs, Scores};

#[derive(Parser, Debug)]
pub struct Args {
    /// A pipeline directory, or the name of one under benchmarks/mgy/.
    #[arg(value_name = "recall|cloud-search|hit-loss")]
    pipeline: String,

    /// The per-family cutoffs a calibration learned. Defaults to the committed
    /// ones.
    #[arg(long, value_name = "cutoffs.txt")]
    cutoffs: Option<PathBuf>,

    /// Which decoy to cut at. The cutoffs file holds each family's five
    /// best-scoring decoys, so this admits at most `c` false positives per
    /// family.
    #[arg(short = 'c', default_value_t = 2, value_name = "N")]
    c: usize,

    /// The query set that was searched. Defaults to the shared one.
    #[arg(long, value_name = "query.hmm")]
    queries: Option<PathBuf>,

    /// The directory holding the target shards. Defaults to the shared one.
    #[arg(long, value_name = "dir")]
    targets: Option<PathBuf>,

    #[arg(short, long, value_name = "scores.tbl")]
    out: Option<PathBuf>,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let dir = match PathBuf::from(&args.pipeline) {
        path if path.is_dir() => path,
        _ => crate::dir().join(&args.pipeline),
    };

    if !dir.is_dir() {
        bail!("no pipeline directory at {}", dir.display());
    }

    let cutoffs = match args.cutoffs {
        Some(path) => path,
        None => tools::mgy_cutoffs()?,
    };

    let query_hmm = args
        .queries
        .unwrap_or_else(|| crate::queries().join("query.hmm"));
    let targets = args.targets.unwrap_or_else(crate::targets);

    let out = args.out.unwrap_or_else(|| dir.join("scores.tbl"));

    let scores = Scores::collect(
        &dir,
        Inputs {
            query_hmm: &query_hmm,
            targets: &targets,
        },
        &cutoffs,
        args.c,
    )?;

    scores.write(&out)?;

    println!(
        "wrote {} ({} pairs across {} runs)",
        out.display(),
        scores.rows.len(),
        scores.runs.len()
    );

    Ok(())
}
