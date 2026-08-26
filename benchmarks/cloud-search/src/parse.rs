//! Turning a finished sweep into numbers, in two steps.
//!
//! `scores` reads every results table once and writes scores.tbl, which is
//! every score everything gave every pair. `grid` reads that back and works out
//! the surface. The split is there because the first half is the expensive one
//! and the one least likely to change: a different statistic, or a different
//! shape for grid.tbl, is a re-parse rather than a re-run of the benchmark.

use std::path::{Path, PathBuf};

use anyhow::Context;
use clap::{Parser, Subcommand};

use crate::cell::Cell;
use crate::scores::{Row, Scores};

#[derive(Subcommand)]
pub enum Cmd {
    /// Read every results table into scores.tbl, one row per pair.
    Scores(ScoresArgs),
    /// Work the surface out of scores.tbl.
    Grid(GridArgs),
}

pub fn main(cmd: Cmd) -> anyhow::Result<()> {
    match cmd {
        Cmd::Scores(args) => scores(args),
        Cmd::Grid(args) => grid(args),
    }
}

#[derive(Parser, Debug)]
pub struct ScoresArgs {
    /// A benchmark directory `cloud-search run` has finished on.
    bench_dir: PathBuf,

    /// The per-family cutoffs mgnify learned. Defaults to its committed ones.
    #[arg(long, value_name = "cutoffs.txt")]
    cutoffs: Option<PathBuf>,

    /// Which decoy to cut at. The cutoffs file holds each family's five
    /// best-scoring decoys, so this admits at most `c` false positives per
    /// family. It is fixed here rather than in `grid`, since the cutoff travels
    /// in the table.
    #[arg(short = 'c', default_value_t = 2, value_name = "N")]
    c: usize,

    #[arg(short, long, value_name = "scores.tbl")]
    out: Option<PathBuf>,
}

fn scores(args: ScoresArgs) -> anyhow::Result<()> {
    let cutoffs = match args.cutoffs {
        Some(path) => path,
        None => tools::mgy_cutoffs()?,
    };
    let out = args
        .out
        .unwrap_or_else(|| args.bench_dir.join("scores.tbl"));

    let scores = Scores::collect(&args.bench_dir, &cutoffs, args.c)?;
    scores.write(&out)?;

    println!(
        "wrote {} ({} pairs across {} cells)",
        out.display(),
        scores.rows.len(),
        scores.columns.len()
    );

    Ok(())
}

#[derive(Parser, Debug)]
pub struct GridArgs {
    /// The scores.tbl `parse scores` wrote, or the benchmark directory holding
    /// one.
    scores: PathBuf,

    #[arg(short, long, value_name = "grid.tbl")]
    out: Option<PathBuf>,
}

fn grid(args: GridArgs) -> anyhow::Result<()> {
    let path = match args.scores.is_dir() {
        true => args.scores.join("scores.tbl"),
        false => args.scores.clone(),
    };

    let scores = Scores::read(&path)?;

    let out = args
        .out
        .unwrap_or_else(|| path.parent().unwrap_or(Path::new(".")).join("grid.tbl"));

    // what hmmer found, at the same cutoff every cell is held to. this is the
    // denominator, so it is worked out the same way a cell's hits are.
    let hmmer_hits = scores.rows.iter().filter(|r| r.clears(r.hmmer())).count();

    anyhow::ensure!(
        hmmer_hits > 0,
        "hmmer found nothing that clears a cutoff; there is nothing to measure against"
    );

    let rows: Vec<GridRow> = scores
        .columns
        .iter()
        .enumerate()
        .map(|(i, column)| {
            let hits = scores
                .rows
                .iter()
                .filter(|r| found(r, i) && r.clears(r.hmmer()))
                .count();

            GridRow {
                cell: column.cell,
                wall_s: column.wall_s,
                hits,
                sens: hits as f64 / hmmer_hits as f64,
            }
        })
        .collect();

    write(&out, &rows, &scores, hmmer_hits)?;

    println!("wrote {}", out.display());
    Ok(())
}

/// Whether one cell found a pair and scored it over the cutoff.
fn found(row: &Row, column: usize) -> bool {
    row.clears(row.scores[column])
}

struct GridRow {
    cell: Cell,
    wall_s: f64,
    hits: usize,
    sens: f64,
}

fn write(path: &Path, rows: &[GridRow], scores: &Scores, hmmer_hits: usize) -> anyhow::Result<()> {
    let headers = ["A", "B", "wall_s", "hits", "sens"];

    let cells: Vec<[String; 5]> = rows
        .iter()
        .map(|r| {
            let (a, b) = match r.cell {
                Cell::Pruned { a, b } => (format!("{a:.1}"), format!("{b:.1}")),
                Cell::Full => ("-".to_string(), "-".to_string()),
            };

            [
                a,
                b,
                format!("{:.4}", r.wall_s),
                r.hits.to_string(),
                format!("{:.4}", r.sens),
            ]
        })
        .collect();

    let widths: Vec<usize> = headers
        .iter()
        .enumerate()
        .map(|(i, h)| {
            cells
                .iter()
                .map(|c| c[i].len())
                .chain(std::iter::once(h.len()))
                .max()
                .unwrap_or(0)
        })
        .collect();

    let mut out = String::new();

    // what was searched, then what the sens column is a fraction of and the two
    // times the figures want as reference lines
    out.push_str(&format!(
        "# query  {:>9} families  {:>12} residues  {:>12} bytes\n\
         # target {:>9} seqs      {:>12} residues  {:>12} bytes\n\
         # pairs  {:>9} rows      {:>12} cells\n\
         # hmmer  {:>9} hits      {:>12.4} wall_s\n\
         # seed   {:>9}           {:>12.4} wall_s\n\
         #\n",
        scores.query.count,
        scores.query.residues,
        scores.query.bytes,
        scores.target.count,
        scores.target.residues,
        scores.target.bytes,
        scores.rows.len(),
        scores.columns.len(),
        hmmer_hits,
        scores.hmmer_wall_s,
        "",
        scores.seed_wall_s,
    ));

    out.push('#');
    for (h, &w) in headers.iter().zip(&widths) {
        out.push_str(&format!(" {h:<w$}"));
    }

    out.push_str("\n#");
    for &w in &widths {
        out.push_str(&format!(" {}", "-".repeat(w)));
    }
    out.push('\n');

    for row in &cells {
        // the two the `# ` takes on a header line, so the columns sit under
        // their names rather than beside them
        out.push_str("  ");
        for (i, (c, &w)) in row.iter().zip(&widths).enumerate() {
            if i > 0 {
                out.push(' ');
            }
            out.push_str(&format!("{c:<w$}"));
        }
        out.push('\n');
    }

    std::fs::write(path, out).with_context(|| format!("failed to write {}", path.display()))
}
