//! Turning any finished pipeline into the one table, and that table into the
//! numbers.
//!
//! There is nothing per-pipeline in here. What the columns are comes out of
//! `manifest.tbl`, so `parse scores recall` and `parse scores cloud-search`
//! are the same code pointed at different directories -- and the analyses that
//! follow are the same code again, pointed at what it wrote.

use std::path::{Path, PathBuf};

use anyhow::{Context, bail};
use clap::{Parser, Subcommand};

use crate::analyze;
use crate::inputs;
use crate::scores::{Inputs, Scores};

#[derive(Subcommand)]
pub enum Cmd {
    /// Read every results table into scores.tbl, one row per pair.
    Scores(ScoresArgs),
    /// What every run found and what it cost, one row per run.
    Summary(TableArgs),
    /// Where the hits hmmer found were lost, one row per checkpoint.
    Funnel(TableArgs),
}

pub fn main(cmd: Cmd) -> anyhow::Result<()> {
    match cmd {
        Cmd::Scores(args) => scores(args),
        Cmd::Summary(args) => derive(args, "summary.tbl", analyze::summary),
        Cmd::Funnel(args) => derive(args, "funnel.tbl", analyze::funnel),
    }
}

#[derive(Parser, Debug)]
pub struct ScoresArgs {
    /// A pipeline directory, or the name of one under benchmarks/mgy/runs/.
    #[arg(value_name = "recall|cloud-search|hit-loss")]
    pipeline: String,

    /// The per-family cutoffs a calibration learned. Defaults to the committed
    /// ones.
    #[arg(long, value_name = "cutoffs.txt")]
    cutoffs: Option<PathBuf>,

    /// Which decoy to cut at. The cutoffs file holds each family's five
    /// best-scoring decoys, so this admits at most `c` false positives per
    /// family. It is fixed here rather than in the analyses, since the cutoff
    /// travels in the table.
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

fn scores(args: ScoresArgs) -> anyhow::Result<()> {
    let dir = pipeline(&args.pipeline)?;

    let cutoffs = match args.cutoffs {
        Some(path) => path,
        None => tools::mgy_cutoffs()?,
    };

    let query_hmm = args.queries.unwrap_or_else(inputs::fixed::query_hmm);
    let targets = args.targets.unwrap_or_else(inputs::fixed::targets);

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

#[derive(Parser, Debug)]
pub struct TableArgs {
    /// The scores.tbl `parse scores` wrote, the pipeline directory holding
    /// one, or the name of one under benchmarks/mgy/runs/.
    #[arg(value_name = "recall|cloud-search|hit-loss")]
    scores: String,

    #[arg(short, long, value_name = "out.tbl")]
    out: Option<PathBuf>,
}

/// Reads scores.tbl back and hands it to one of the analyses.
fn derive(
    args: TableArgs,
    name: &str,
    f: fn(&Scores, &Path) -> anyhow::Result<()>,
) -> anyhow::Result<()> {
    let path = match PathBuf::from(&args.scores) {
        path if path.is_file() => path,
        _ => pipeline(&args.scores)?.join("scores.tbl"),
    };

    if !path.is_file() {
        bail!(
            "no scores.tbl at {}; run `mgy parse scores` first",
            path.display()
        );
    }

    let out = match args.out {
        Some(path) => path,
        None => path
            .parent()
            .context("scores.tbl has no directory")?
            .join(name),
    };

    f(&Scores::read(&path)?, &out)?;

    println!("wrote {}", out.display());
    Ok(())
}

/// A pipeline directory, given either as a path or by name.
fn pipeline(name: &str) -> anyhow::Result<PathBuf> {
    let dir = match PathBuf::from(name) {
        path if path.is_dir() => path,
        _ => crate::runs().join(name),
    };

    if !dir.is_dir() {
        bail!("no pipeline directory at {}", dir.display());
    }

    Ok(dir)
}
