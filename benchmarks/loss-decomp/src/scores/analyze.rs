//! The analyses, which are groupings over a table of pairs.
//!
//! Nothing in here reads a results file. `stages` counts per checkpoint
//! within a run, held to one denominator: what the ceiling kept over its
//! family's cutoff.
//!
//! It is separate from `parse` because reading every results table is the
//! expensive half and the half least likely to change: a different statistic
//! is a re-run of this, not of the benchmark.

use std::path::Path;

use anyhow::{Context, ensure};

use crate::scores::Run;
use crate::scores::runs::Reader as Runs;

/// What was searched, and what the fractions are fractions of.
fn preamble(
    table: &mut toil::Table,
    meta: &crate::scores::Meta,
    ceiling: &Run,
    hits: usize,
    rows: u64,
) {
    let (mut count, mut residues, mut bytes) = (0usize, 0u64, 0u64);
    for (_, size) in &meta.targets {
        count += size.count;
        residues += size.residues;
        bytes += size.bytes;
    }

    table
        .comment(format!(
            "query   {:>9} families  {:>12} residues  {:>12} bytes",
            meta.query.count, meta.query.residues, meta.query.bytes,
        ))
        .comment(format!(
            "target  {count:>9} seqs      {residues:>12} residues  {bytes:>12} bytes"
        ))
        .comment(format!(
            "pairs   {rows:>9} rows      {:>12} runs",
            meta.runs.len()
        ))
        .comment(format!(
            "ceiling {hits:>9} hits      {:>12.4} wall_s  {}",
            ceiling.wall_s, ceiling.name
        ))
        .comment("");
}

/// What one unit's pairs came to, per run.
struct Tally {
    /// Pairs the ceiling kept at or above their family's cutoff, for this
    /// unit.
    ceiling: usize,
    lost_seed: Vec<usize>,
    lost_cloud_align: Vec<usize>,
    /// Scored, and below the family's cutoff. nail found the pair and would
    /// not return it, which is a loss like any other.
    lost_cutoff: Vec<usize>,
    kept: Vec<usize>,
}

impl Tally {
    fn new(runs: usize) -> Tally {
        Tally {
            ceiling: 0,
            lost_seed: vec![0; runs],
            lost_cloud_align: vec![0; runs],
            lost_cutoff: vec![0; runs],
            kept: vec![0; runs],
        }
    }
}

/// Where the hits the ceiling kept are lost, per unit and per run.
pub fn stages(path: &Path, out: &Path) -> anyhow::Result<()> {
    let mut scores = Runs::open(path)?;

    let ceiling = scores.meta().ceiling()?;
    let runs = scores.meta().runs.len();

    ensure!(
        scores.meta().runs.iter().all(|run| run.seeds.is_some()),
        "a run kept no seed list, so there is no seeding checkpoint to split on"
    );

    // per unit as well as per run. A `cross` set searches one query against
    // several kinds of target, and what the ceiling kept in one is not the
    // denominator for another: summed, two corpora make a sensitivity that
    // describes neither
    let mut at: indexmap::IndexMap<String, Tally> = indexmap::IndexMap::new();
    let mut rows = 0u64;

    scores.each(|row| {
        rows += 1;

        if !row.row().passed(ceiling) {
            return Ok(());
        }

        let unit = row.row().shard().to_string();
        let tally = at.entry(unit).or_insert_with(|| Tally::new(runs));
        tally.ceiling += 1;

        for run in 0..runs {
            if run == ceiling {
                continue;
            }

            // per run rather than per pair: every arm has its own seed list,
            // so whether the pair was ever offered is the arm's answer
            match (row.seeded(run), row.present(run), row.row().passed(run)) {
                (false, _, _) => tally.lost_seed[run] += 1,
                (_, false, _) => tally.lost_cloud_align[run] += 1,
                (_, _, false) => tally.lost_cutoff[run] += 1,
                (_, _, true) => tally.kept[run] += 1,
            }
        }

        Ok(())
    })?;

    ensure!(
        at.values().any(|t| t.ceiling > 0),
        "the ceiling kept nothing over a cutoff; there is nothing to measure against"
    );

    // one row per (unit, run), a column per checkpoint. the ceiling has no
    // row: it loses nothing of its own by definition, and its count is on
    // every row of its unit
    let headers = [
        "unit",
        "run",
        "ceiling",
        "lost_seed",
        "lost_align",
        "lost_cutoff",
        "kept",
        "sens",
    ]
    .map(str::to_string)
    .to_vec();
    let mut cells: Vec<Vec<String>> = Vec::new();

    for (unit, tally) in &at {
        for (run, column) in scores.meta().runs.iter().enumerate() {
            if run == ceiling {
                continue;
            }

            cells.push(vec![
                unit.clone(),
                column.name.clone(),
                tally.ceiling.to_string(),
                tally.lost_seed[run].to_string(),
                tally.lost_cloud_align[run].to_string(),
                tally.lost_cutoff[run].to_string(),
                tally.kept[run].to_string(),
                format!("{:.4}", frac(tally.kept[run], tally.ceiling)),
            ]);
        }
    }

    ensure!(
        !cells.is_empty(),
        "nothing but the ceiling ran, so there is no arm to measure"
    );

    let hits: usize = at.values().map(|tally| tally.ceiling).sum();
    let top = &scores.meta().runs[ceiling];

    let mut table = toil::Table::new(toil::Schema::new(headers));
    preamble(&mut table, scores.meta(), top, hits, rows);
    table.meta("ceiling", [top.name.as_str()]);
    for (unit, tally) in &at {
        table.meta("hits", [unit.as_str(), &tally.ceiling.to_string()]);
    }
    for row in cells {
        table.row(row);
    }

    table
        .write(out)
        .with_context(|| format!("failed to write {}", out.display()))
}

fn frac(n: usize, of: usize) -> f64 {
    match of {
        0 => 0.0,
        of => n as f64 / of as f64,
    }
}
