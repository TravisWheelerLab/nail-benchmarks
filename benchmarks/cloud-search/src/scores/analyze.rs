//! The summary, which is a grouping over a table of pairs.
//!
//! Nothing in here reads a results file. A summary counts per run what that
//! run reported over nail's cutoff, and holds it against what the full run
//! reported. A run is one cell over all its reps: the hits are rep 1's, and
//! the wall clock is the median over the reps, with the fastest and slowest
//! beside it.
//!
//! A family with no nail cutoff is left out of every count. The table says how many families that was in a
//! `#= left_out families <n>` line, and counts per run the pairs that run
//! reported in them in a `left_out` column.
//!
//! It is separate from `parse` because reading every results table is the
//! expensive half and the half least likely to change: a different statistic
//! is a re-run of this, not of the benchmark.

use std::collections::{BTreeSet, HashSet};
use std::path::Path;

use anyhow::{Context, ensure};

use crate::run::FULL;
use crate::scores::complete;
use crate::scores::runs::Reader as Runs;

/// What every run found, and what it cost.
///
/// One pass over the table, counting per run what it reported over the
/// cutoff. Nothing is held but the counters.
pub fn summary(path: &Path, out: &Path) -> anyhow::Result<()> {
    let mut scores = Runs::open(path)?;

    let runs = scores.meta().runs.len();
    let complete = complete(&scores.meta().cutoffs, scores.meta().c)?;

    let full = scores
        .meta()
        .runs
        .iter()
        .position(|run| run.name == FULL)
        .with_context(|| format!("no `{FULL}` run in {} to measure against", path.display()))?;

    let mut hits = vec![0usize; runs];
    let mut left = vec![0usize; runs];
    let mut left_families = HashSet::new();
    let mut rows = 0u64;

    scores.each(|scored| {
        let row = scored.row();

        if left_out(&complete, &mut left_families, row.field(0))? {
            for (run, n) in left.iter_mut().enumerate() {
                *n += usize::from(scored.present(run));
            }
            return Ok(());
        }

        rows += 1;

        for (run, n) in hits.iter_mut().enumerate() {
            *n += usize::from(row.passed(run));
        }

        Ok(())
    })?;

    ensure!(
        hits[full] > 0,
        "the `{FULL}` run found nothing that clears a cutoff; there is nothing to measure against"
    );

    // every setting any run recorded, so a pipeline that swept two knobs gets
    // two columns and one that swept none gets none
    let meta = scores.meta();
    let keys: BTreeSet<&str> = meta
        .runs
        .iter()
        .flat_map(|run| run.params.keys())
        .map(String::as_str)
        .collect();

    let mut headers = vec!["name".to_string()];
    headers.extend(keys.iter().map(|k| k.to_string()));
    headers.extend(
        [
            "reps", "wall_s", "wall_min", "wall_max", "hits", "coverage", "left_out",
        ]
        .iter()
        .map(|h| h.to_string()),
    );

    let mut table = toil::Table::new(toil::Schema::new(headers));
    table.meta(
        "left_out",
        [toil::Cell::from("families"), left_families.len().into()],
    );
    preamble(
        &mut table,
        meta,
        hits[full],
        median(&meta.runs[full].walls),
        rows,
    );

    for (i, run) in meta.runs.iter().enumerate() {
        let mut cells = vec![run.name.clone()];
        cells.extend(keys.iter().map(|k| match run.params.get(*k) {
            Some(value) => value.clone(),
            None => "-".to_string(),
        }));
        let min = run.walls.iter().copied().fold(f64::INFINITY, f64::min);
        let max = run.walls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        cells.extend([
            run.walls.len().to_string(),
            format!("{:.4}", median(&run.walls)),
            format!("{min:.4}"),
            format!("{max:.4}"),
            hits[i].to_string(),
            format!("{:.4}", hits[i] as f64 / hits[full] as f64),
            left[i].to_string(),
        ]);

        table.row(cells);
    }

    table
        .write(out)
        .with_context(|| format!("failed to write {}", out.display()))
}

/// The middle of `walls`, or the mean of the middle two.
fn median(walls: &[f64]) -> f64 {
    let mut sorted = walls.to_vec();
    sorted.sort_by(f64::total_cmp);

    let mid = sorted.len() / 2;
    match sorted.len() % 2 {
        1 => sorted[mid],
        _ => (sorted[mid - 1] + sorted[mid]) / 2.0,
    }
}

/// What was searched, what the fractions are fractions of, and the two times
/// the figures use as reference lines.
fn preamble(
    table: &mut toil::Table,
    meta: &crate::scores::Meta,
    full: usize,
    wall: f64,
    rows: u64,
) {
    let (mut count, mut residues, mut bytes) = (0usize, 0u64, 0u64);
    for (_, size) in &meta.targets {
        count += size.count;
        residues += size.residues;
        bytes += size.bytes;
    }

    // a dash rather than a zero for a pipeline that never seeded: seeding
    // taking no time and there being no seeding are different things
    let seed = match meta.seeds.is_empty() {
        true => "-".to_string(),
        false => format!("{:.4}", meta.seeds.iter().map(|(_, w)| w).sum::<f64>()),
    };

    table
        .comment(format!(
            "query  {:>9} families  {:>12} residues  {:>12} bytes",
            meta.query.count, meta.query.residues, meta.query.bytes,
        ))
        .comment(format!(
            "target {count:>9} seqs      {residues:>12} residues  {bytes:>12} bytes"
        ))
        .comment(format!(
            "pairs  {rows:>9} rows      {:>12} runs",
            meta.runs.len()
        ))
        .comment(format!("full   {full:>9} hits      {wall:>12.4} wall_s"))
        .comment(format!("seed   {:>9}           {seed:>12} wall_s", ""))
        .comment("");
}

/// Whether a row's family is left out, adding it to `met` when it is.
fn left_out(
    complete: &HashSet<String>,
    met: &mut HashSet<String>,
    query: &[u8],
) -> anyhow::Result<bool> {
    let query = std::str::from_utf8(query).context("a query name that is not utf-8")?;

    if complete.contains(query) {
        return Ok(false);
    }

    if !met.contains(query) {
        met.insert(query.to_string());
    }

    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// What `analysis` writes from the fixture table.
    fn written(name: &str, analysis: fn(&Path, &Path) -> anyhow::Result<()>) -> String {
        written_from(name, &crate::scores::runs::tests::file(), analysis)
    }

    /// What `analysis` writes from `table`.
    fn written_from(
        name: &str,
        table: &str,
        analysis: fn(&Path, &Path) -> anyhow::Result<()>,
    ) -> String {
        let dir = std::env::temp_dir().join(format!("mgy-analyze-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let path = dir.join(format!("{name}.runs.tbl"));
        std::fs::write(&path, table).unwrap();

        let out = dir.join(format!("{name}.tbl"));
        analysis(&path, &out).unwrap();
        std::fs::read_to_string(&out).unwrap()
    }

    /// The last cell of the first row that opens with the cells `key`.
    fn left_column(text: &str, key: &[&str]) -> String {
        text.lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split_whitespace().collect::<Vec<_>>())
            .find(|cells| cells.starts_with(key))
            .and_then(|cells| cells.last().map(|cell| cell.to_string()))
            .unwrap()
    }

    #[test]
    fn a_summary_leaves_out_a_family_short_of_a_cutoff_and_counts_it() {
        let text = written("summary", summary);

        assert!(text.contains("#= left_out families 1\n"));

        assert_eq!(left_column(&text, &["A2.0-B4.0"]), "1");
        assert_eq!(left_column(&text, &["full"]), "1");
        assert!(text.contains("pairs          3 rows"));
    }

    #[test]
    fn coverage_is_a_run_s_hits_over_the_full_run_s() {
        let text = written("coverage", summary);

        // the full run cleared alpha, beta's first pair and beta's second;
        // the cell cleared only beta's second
        assert!(text.contains("full           3 hits"));

        let cells: Vec<&str> = text
            .lines()
            .find(|line| line.starts_with("A2.0-B4.0"))
            .unwrap()
            .split_whitespace()
            .collect();
        assert_eq!(&cells[cells.len() - 3..cells.len() - 1], ["1", "0.3333"]);
    }

    /// The cells of the first row that opens with `name`.
    fn row<'t>(text: &'t str, name: &str) -> Vec<&'t str> {
        text.lines()
            .find(|line| line.split_whitespace().next() == Some(name))
            .unwrap()
            .split_whitespace()
            .collect()
    }

    #[test]
    fn the_wall_clock_is_the_median_over_reps() {
        let text = written("median", summary);

        // the fixture's cell ran three reps at 3, 5 and 4 seconds
        let cells = row(&text, "A2.0-B4.0");
        assert_eq!(&cells[3..7], ["3", "4.0000", "3.0000", "5.0000"]);

        // and the full run two, at 9 and 10
        let cells = row(&text, "full");
        assert_eq!(&cells[3..7], ["2", "9.5000", "9.0000", "10.0000"]);
        let preamble = text.lines().find(|l| l.starts_with("# full")).unwrap();
        assert!(preamble.ends_with(" 9.5000 wall_s"), "{preamble}");
    }

    #[test]
    fn a_table_of_one_rep_still_summarizes() {
        let one = crate::scores::runs::tests::file()
            .replace("3.0000,5.0000,4.0000", "3.0000")
            .replace("9.0000,10.0000", "9.0000");
        let text = written_from("one-rep", &one, summary);

        let cells = row(&text, "A2.0-B4.0");
        assert_eq!(&cells[3..7], ["1", "3.0000", "3.0000", "3.0000"]);
    }
}
