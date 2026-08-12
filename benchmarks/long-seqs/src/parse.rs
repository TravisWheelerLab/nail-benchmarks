use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::Context;
use clap::{Parser, Subcommand};

use crate::{dir, PAIRS};

/// Analysis subcommands for this benchmark.
#[derive(Subcommand)]
pub enum Cmd {
    /// Emit `cells.long.txt`: DP matrix area against the fraction of it
    /// computed, one row per pair. The pct-id benchmark overlays this on its
    /// own cell-fraction figure via `--long_hits`.
    Cells(CellsArgs),
}

#[derive(Parser)]
pub struct CellsArgs {
    /// Results directory; defaults to this benchmark's results/.
    #[arg(short, long)]
    results: Option<PathBuf>,

    /// Where to write cells.long.txt.
    #[arg(short, long)]
    out: Option<PathBuf>,

    /// Name of the run to read, matching bench.toml.
    #[arg(long, default_value = "nail.seq")]
    run: String,

}

pub fn main(cmd: Cmd) -> anyhow::Result<()> {
    match cmd {
        Cmd::Cells(args) => cells(args),
    }
}

fn cells(args: CellsArgs) -> anyhow::Result<()> {
    let dir = dir();

    let results = args.results.unwrap_or_else(|| dir.join("results"));
    let out_path = args.out.unwrap_or_else(|| results.join("cells.long.txt"));

    let mut out = BufWriter::new(
        File::create(&out_path)
            .with_context(|| format!("failed to create {}", out_path.display()))?,
    );

    // resolve output paths through the runs table rather than rebuilding the
    // naming convention here
    let runs = run::Runs::from_dir(&results)?;

    for i in 1..=PAIRS {
        let tbl = runs.table_path(&args.run, &i.to_string());
        let cell_frac = last_cell_frac(&tbl)
            .with_context(|| format!("failed to read a hit from {}", tbl.display()))?;

        let q_len = bioio::fasta::residue_len(dir.join(format!("query/{i}.query.fa")))?;
        let t_len = bioio::fasta::residue_len(dir.join(format!("target/{i}.target.fa")))?;

        writeln!(out, "{},{:.5}", q_len * t_len, cell_frac)?;
    }

    out.flush()?;
    println!("wrote {}", out_path.display());
    Ok(())
}

/// Cell fraction of the last hit in a nail table, which is the one these
/// single-pair searches are about.
fn last_cell_frac(path: &Path) -> anyhow::Result<f64> {
    let tbl = bioio::tbl::nail::NailTable::from_path(path)
        .with_context(|| format!("failed to read {}", path.display()))?;

    tbl.hits
        .last()
        .map(|h| h.cell_frac)
        .context("no hits in table")
}

