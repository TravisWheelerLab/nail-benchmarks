//! Turning a finished run into the file the cell-fraction figure overlays.
//!
//! What ran comes out of `ledger.tbl` -- which pair each search covered and
//! where it put its table -- rather than out of a count of the pairs that are
//! checked in. A pair whose search failed is left out with a warning instead of
//! coming back as a missing file.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};
use clap::{Parser, Subcommand};

use util::ledger::{self, Ledger};
use util::manifest;

use util::set::Set;

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
    /// Where the run wrote. Defaults to this benchmark's outputs/
    #[arg(long, value_name = "dir")]
    out: Option<PathBuf>,

    /// Where cells.long.txt goes. Defaults to figures/ beside the run
    #[arg(short, long, value_name = "dir")]
    figures: Option<PathBuf>,

    /// Which run's tables to read cell fractions from
    #[arg(long, value_name = "NAME", default_value = RUN_NAME)]
    run: String,

    /// Which label of paths.toml to read. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,
}

pub fn main(cmd: Cmd, paths: &crate::Paths) -> anyhow::Result<()> {
    match cmd {
        Cmd::Cells(args) => cells(args, paths),
    }
}

fn cells(args: CellsArgs, paths: &crate::Paths) -> anyhow::Result<()> {
    let out = args.out.unwrap_or_else(|| paths.run.clone());
    let figures = args.figures.unwrap_or_else(|| paths.analysis.clone());

    let set = Set::load_as(&paths.set, &util::set::shape::PAIRS)?;

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

        // the lengths the set wrote down, rather than the files read again
        let unit = set
            .units()
            .find(|u| u.name() == pair)
            .with_context(|| format!("the set has no unit {pair:?}"))?;

        let q_len = unit.number("query_residues")?;
        let t_len = unit.number("residues")?;

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
/// error rather than something to skip: it would mean the record and this
/// analysis disagree about what was measured.
fn searches(out: &Path, run: &str) -> anyhow::Result<Vec<(String, PathBuf)>> {
    let ran = Ledger::load(out)?;
    ledger::warn(ran.failed(), "pair(s)");

    let results = out.join("results");
    let mut searches = Vec::new();

    for row in ran.runs().filter(|row| row.name == run) {
        ensure!(
            row.tool == "nail",
            "run {run:?} was produced by {}, which reports no cell fractions",
            row.tool
        );
        ensure!(
            !row.shard.is_empty(),
            "run {run:?} has no shard saying which pair it was"
        );

        searches.push((
            row.shard.clone(),
            manifest::table_path(&results, run, &row.shard),
        ));
    }

    if searches.is_empty() {
        bail!("no finished {run:?} runs in {}", out.display());
    }

    Ok(searches)
}

/// Cell fraction of the last hit in a nail table, which is the one these
/// single-pair searches are about.
fn last_cell_frac(path: &Path) -> anyhow::Result<f64> {
    util::nail::cell_fracs(path)?
        .last()
        .map(|h| h.cell_frac)
        .context("no hits in table")
}
