//! Turning a finished run into the file the cell-fraction figure overlays.
//!
//! What ran comes out of `manifest.tbl` -- which pair each search covered and
//! where it put its table -- rather than out of a count of the pairs that are
//! checked in. A pair whose search failed is left out with a warning instead of
//! coming back as a missing file.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};
use clap::{Parser, Subcommand};
use libsail::collection::Indexable;
use libsail::seq::fasta::Fasta;

use bench::manifest::{self, Manifest};

use crate::inputs;
use crate::run::RUN_NAME;

/// Analysis subcommands for this benchmark.
#[derive(Subcommand)]
pub enum Cmd {
    /// Emit `cells.long.txt`: DP matrix area against the fraction of it
    /// computed, one row per pair. The pid benchmark overlays this on its
    /// own cell-fraction figure via `--long_hits`.
    Cells(CellsArgs),
}

#[derive(Parser)]
pub struct CellsArgs {
    /// Where the run wrote. Defaults to this benchmark's outputs/.
    #[arg(long, value_name = "dir")]
    out: Option<PathBuf>,

    /// Where cells.long.txt goes. Defaults to figures/ beside the run.
    #[arg(short, long, value_name = "dir")]
    figures: Option<PathBuf>,

    /// Which run's tables to read cell fractions from.
    #[arg(long, value_name = "NAME", default_value = RUN_NAME)]
    run: String,
}

pub fn main(cmd: Cmd) -> anyhow::Result<()> {
    match cmd {
        Cmd::Cells(args) => cells(args),
    }
}

fn cells(args: CellsArgs) -> anyhow::Result<()> {
    let out = args.out.unwrap_or_else(inputs::outputs);
    let figures = args.figures.unwrap_or_else(|| out.join("figures"));

    let searches = searches(&out, &args.run)?;

    std::fs::create_dir_all(&figures)
        .with_context(|| format!("failed to create {}", figures.display()))?;

    let path = figures.join("cells.long.txt");
    let mut file = BufWriter::new(
        File::create(&path).with_context(|| format!("failed to create {}", path.display()))?,
    );

    for (pair, table) in searches {
        let cell_frac = last_cell_frac(&table)
            .with_context(|| format!("failed to read a hit from {}", table.display()))?;

        let q_len = residue_len(&inputs::query(&pair))?;
        let t_len = residue_len(&inputs::target(&pair))?;

        writeln!(file, "{},{:.5}", q_len * t_len, cell_frac)?;
    }

    file.flush()?;
    println!("wrote {}", path.display());
    Ok(())
}

/// The pairs one run covered and the table it wrote for each, in the order the
/// pipeline declared them.
///
/// Only nail reports a cell fraction, so a row filed under another tool is an
/// error rather than something to skip: it would mean the manifest and this
/// analysis disagree about what was measured.
fn searches(out: &Path, run: &str) -> anyhow::Result<Vec<(String, PathBuf)>> {
    let manifest = Manifest::read(&out.join("manifest.tbl"))?;
    let results = out.join("results");

    let failed: Vec<&str> = manifest
        .failed()
        .filter_map(|row| row.get(manifest::SHARD))
        .collect();
    if !failed.is_empty() {
        eprintln!(
            "warning: leaving out {} pair(s) that did not finish: {}",
            failed.len(),
            failed.join(", ")
        );
    }

    let mut searches = Vec::new();

    for row in manifest.runs() {
        if row.get(manifest::NAME) != Some(run) {
            continue;
        }

        let tool = row
            .get(manifest::TOOL)
            .with_context(|| format!("run {run:?} has no tool"))?;
        ensure!(
            tool == "nail",
            "run {run:?} was produced by {tool}, which reports no cell fractions"
        );

        let pair = row
            .get(manifest::SHARD)
            .with_context(|| format!("run {run:?} has no shard saying which pair it was"))?;

        searches.push((pair.to_string(), manifest::table_path(&results, run, pair)));
    }

    if searches.is_empty() {
        bail!("no finished {run:?} rows in {}/manifest.tbl", out.display());
    }

    Ok(searches)
}

/// Cell fraction of the last hit in a nail table, which is the one these
/// single-pair searches are about.
fn last_cell_frac(path: &Path) -> anyhow::Result<f64> {
    bench::nail::cell_fracs(path)?
        .last()
        .map(|h| h.cell_frac)
        .context("no hits in table")
}

/// Residue count of a fasta holding exactly one sequence.
fn residue_len(path: &Path) -> anyhow::Result<usize> {
    let fa =
        Fasta::open(path).with_context(|| format!("failed to parse {}", path.display()))?;

    ensure!(
        fa.len() == 1,
        "expected exactly one sequence in {}, found {}",
        path.display(),
        fa.len()
    );

    Ok(fa.get(0).expect("one record").seq.len())
}
