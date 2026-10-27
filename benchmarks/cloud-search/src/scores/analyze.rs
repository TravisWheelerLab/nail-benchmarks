//! The analyses, which are groupings over a table of pairs.
//!
//! Nothing in here reads a results file. A summary counts per run, and stages
//! counts per checkpoint within a run; both are held to the same denominator,
//! which is what hmmer found and scored over its family's cutoff.
//!
//! A family short of a cutoff for any tool, or for hmmer's best domain, is
//! left out of every count. Each table says how many families that was in a
//! `#= left_out families <n>` line, and counts per run the pairs that run
//! reported in them in a `left_out` column.
//!
//! They are separate from `parse` because reading every results table is the
//! expensive half and the half least likely to change: a different statistic
//! is a re-run of this, not of the benchmark.

use std::collections::{BTreeSet, HashSet};
use std::path::Path;

use anyhow::{Context, ensure};

use crate::scores::runs::Reader as Runs;
use crate::scores::{Tool, complete};

/// What every run found, and what it cost.
///
/// One pass over the table, counting per run: what it reported over the
/// cutoff, how much of that hmmer also found, and how much of that hmmer read
/// as one domain rather than several. Nothing is held but the counters.
///
pub fn summary(path: &Path, out: &Path) -> anyhow::Result<()> {
    let mut scores = Runs::open(path)?;

    let hmmer = scores.meta().hmmer()?;
    let runs = scores.meta().runs.len();
    let complete = complete(&scores.meta().cutoffs, scores.meta().c)?;

    let mut found = vec![0usize; runs];
    let mut hits = vec![0usize; runs];
    let mut hits_sd = vec![0usize; runs];
    let mut left = vec![0usize; runs];
    let mut left_families = HashSet::new();
    let (mut truth, mut truth_sd) = (0usize, 0usize);
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

        let true_hit = row.passed(hmmer);

        // a hit hmmer breaks into one region is a different question from one
        // it breaks into several: the tools disagree most about the second
        let single = true_hit && row.domain_count() == 1;

        if true_hit {
            truth += 1;
            truth_sd += usize::from(single);
        }

        for run in 0..runs {
            if !row.passed(run) {
                continue;
            }

            found[run] += 1;

            // held to hmmer as well, so a run is credited for what it agreed
            // with rather than for everything it scored highly
            if true_hit {
                hits[run] += 1;
                hits_sd[run] += usize::from(single);
            }
        }

        Ok(())
    })?;

    ensure!(
        truth > 0,
        "hmmer found nothing that clears a cutoff; there is nothing to measure against"
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

    let mut headers = vec!["name".to_string(), "tool".to_string()];
    headers.extend(keys.iter().map(|k| k.to_string()));
    headers.extend(
        [
            "wall_s", "found", "hits", "sens", "hits_sd", "sens_sd", "left_out",
        ]
        .iter()
        .map(|h| h.to_string()),
    );

    let cells: Vec<Vec<String>> = meta
        .runs
        .iter()
        .enumerate()
        .map(|(i, run)| {
            let mut cells = vec![run.name.clone(), run.tool.to_string()];
            cells.extend(keys.iter().map(|k| match run.params.get(*k) {
                Some(value) => value.clone(),
                None => "-".to_string(),
            }));
            cells.extend([
                format!("{:.4}", run.wall_s),
                found[i].to_string(),
                hits[i].to_string(),
                format!("{:.4}", frac(hits[i], truth)),
                hits_sd[i].to_string(),
                format!("{:.4}", frac(hits_sd[i], truth_sd)),
                left[i].to_string(),
            ]);

            cells
        })
        .collect();

    let mut table = toil::Table::new(toil::Schema::new(headers));
    table.meta(
        "left_out",
        [toil::Cell::from("families"), left_families.len().into()],
    );
    preamble(&mut table, meta, truth, rows);
    for row in cells {
        table.row(row);
    }

    table
        .write(out)
        .with_context(|| format!("failed to write {}", out.display()))
}

/// What was searched, what the fractions are fractions of, and the two times
/// the figures use as reference lines.
fn preamble(table: &mut toil::Table, meta: &crate::scores::Meta, truth: usize, rows: u64) {
    let (mut count, mut residues, mut bytes) = (0usize, 0u64, 0u64);
    for (_, size) in &meta.targets {
        count += size.count;
        residues += size.residues;
        bytes += size.bytes;
    }

    let hmmer = meta
        .runs
        .iter()
        .find(|run| run.tool == Tool::Hmmer)
        .map(|run| run.wall_s)
        .unwrap_or_default();

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
        .comment(format!("hmmer  {truth:>9} hits      {hmmer:>12.4} wall_s"))
        .comment(format!("seed   {:>9}           {seed:>12} wall_s", ""))
        .comment("");
}

/// What one unit's pairs came to, per run.
struct Tally {
    /// Pairs hmmer reported at or above their family's cutoff, for this unit.
    truth: usize,
    lost_seed: Vec<usize>,
    lost_cloud_align: Vec<usize>,
    /// Scored, and below the family's cutoff. nail found the pair and would
    /// not return it, which is a loss like any other.
    lost_cutoff: Vec<usize>,
    kept: Vec<usize>,
    /// Pairs each run reported in a family left out, for this unit.
    left: Vec<usize>,
}

impl Tally {
    fn new(runs: usize) -> Tally {
        Tally {
            truth: 0,
            lost_seed: vec![0; runs],
            lost_cloud_align: vec![0; runs],
            lost_cutoff: vec![0; runs],
            kept: vec![0; runs],
            left: vec![0; runs],
        }
    }
}

/// Where the hits hmmer found are lost, per unit and per run.
pub fn stages(path: &Path, out: &Path) -> anyhow::Result<()> {
    let mut scores = Runs::open(path)?;

    let hmmer = scores.meta().hmmer()?;
    let runs = scores.meta().runs.len();
    let complete = complete(&scores.meta().cutoffs, scores.meta().c)?;
    let mut left_families = HashSet::new();

    ensure!(
        !scores.meta().seeds.is_empty(),
        "this pipeline kept no seeds, so there is no seeding checkpoint to split on"
    );

    // per unit as well as per run. A `cross` set searches one query against
    // several kinds of target, and what hmmer found in one is not the truth
    // set for another: summed, the two corpora make a sensitivity that
    // describes neither
    let mut at: indexmap::IndexMap<String, Tally> = indexmap::IndexMap::new();
    let mut rows = 0u64;

    scores.each(|row| {
        let unit = row.row().shard();

        if left_out(&complete, &mut left_families, row.row().field(0))? {
            let tally = tally(&mut at, unit, runs);
            for (run, n) in tally.left.iter_mut().enumerate() {
                *n += usize::from(row.present(run));
            }
            return Ok(());
        }

        rows += 1;

        if !row.row().passed(hmmer) {
            return Ok(());
        }

        let tally = tally(&mut at, unit, runs);
        tally.truth += 1;

        for run in 0..runs {
            if run == hmmer {
                continue;
            }

            // per run rather than per pair: a seeding sweep gives every arm
            // its own seed list, so whether the pair was ever offered is the
            // arm's answer and not the pipeline's
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
        at.values().any(|t| t.truth > 0),
        "hmmer found nothing that clears a cutoff; there is nothing to measure against"
    );

    // one row per (unit, run), a column per checkpoint. Written long it was
    // four rows apiece, where `n` meant a population on two of them and a loss
    // on the other two, and the last row's fraction only ever repeated the one
    // above it
    let headers = [
        "unit",
        "run",
        "truth",
        "lost_seed",
        "lost_align",
        "lost_cutoff",
        "kept",
        "sens",
        "left_out",
    ]
    .map(str::to_string)
    .to_vec();
    let mut cells: Vec<Vec<String>> = Vec::new();

    for (unit, tally) in &at {
        for (run, column) in scores.meta().runs.iter().enumerate() {
            if run == hmmer {
                continue;
            }

            cells.push(vec![
                unit.clone(),
                column.name.clone(),
                tally.truth.to_string(),
                tally.lost_seed[run].to_string(),
                tally.lost_cloud_align[run].to_string(),
                tally.lost_cutoff[run].to_string(),
                tally.kept[run].to_string(),
                format!("{:.4}", frac(tally.kept[run], tally.truth)),
                tally.left[run].to_string(),
            ]);
        }
    }

    ensure!(
        !cells.is_empty(),
        "nothing but hmmer ran, so there is no pipeline to trace"
    );

    let truth: usize = at.values().map(|tally| tally.truth).sum();

    let mut table = toil::Table::new(toil::Schema::new(headers));
    table.meta(
        "left_out",
        [toil::Cell::from("families"), left_families.len().into()],
    );
    preamble(&mut table, scores.meta(), truth, rows);
    for row in cells {
        table.row(row);
    }

    table
        .write(out)
        .with_context(|| format!("failed to write {}", out.display()))
}

/// One unit's tally, made on the first row that names the unit.
fn tally<'a>(
    at: &'a mut indexmap::IndexMap<String, Tally>,
    unit: &str,
    runs: usize,
) -> &'a mut Tally {
    if !at.contains_key(unit) {
        at.insert(unit.to_string(), Tally::new(runs));
    }

    &mut at[unit]
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

fn frac(n: usize, of: usize) -> f64 {
    match of {
        0 => 0.0,
        of => n as f64 / of as f64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// What `analysis` writes from the fixture table.
    fn written(name: &str, analysis: fn(&Path, &Path) -> anyhow::Result<()>) -> String {
        let dir = std::env::temp_dir().join(format!("mgy-analyze-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let path = dir.join(format!("{name}.runs.tbl"));
        std::fs::write(&path, crate::scores::runs::tests::file()).unwrap();

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

        // zeta's pair is hmmer's and the first cell's, and not the full
        // run's, which was never offered it
        assert_eq!(left_column(&text, &["hmmer"]), "1");
        assert_eq!(left_column(&text, &["A2.0-B4.0"]), "1");
        assert_eq!(left_column(&text, &["full"]), "0");

        // beta and gamma, and not zeta, which hmmer scored over nail's
        // cutoff
        assert!(text.contains("hmmer          2 hits"));
        assert!(text.contains("pairs          3 rows"));
    }

    #[test]
    fn stages_leave_out_a_family_short_of_a_cutoff_and_count_it() {
        let text = written("stages", stages);

        assert!(text.contains("#= left_out families 1\n"));
        assert_eq!(left_column(&text, &["1", "A2.0-B4.0"]), "1");
        assert_eq!(left_column(&text, &["1", "full"]), "0");
        assert!(text.contains("hmmer          2 hits"));
    }
}
