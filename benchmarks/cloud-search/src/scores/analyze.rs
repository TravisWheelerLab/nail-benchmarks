//! The analyses, which are groupings over a table of pairs.
//!
//! Nothing in here reads a results file. A summary counts per run, and stages
//! counts per checkpoint within a run; both are held to the same denominator,
//! which is what hmmer found and scored over its family's cutoff.
//!
//! They are separate from `parse` because reading every results table is the
//! expensive half and the half least likely to change: a different statistic
//! is a re-run of this, not of the benchmark.

use std::collections::BTreeSet;
use std::path::Path;

use anyhow::{Context, ensure};

use crate::scores::frame::{Frame, Verdict};
use crate::scores::runs::Reader as Runs;
use crate::scores::{Named, Tool, runs};

/// What every run found, and what it cost.
///
/// One pass over the table, counting per run: what it reported over the
/// cutoff, how much of that hmmer also found, and how much of that hmmer read
/// as one domain rather than several. Nothing is held but the counters.
///
/// This reads either table. A summary asks only whether each run cleared its
/// family's cutoff and what the domain list says, and the two tables answer
/// the first differently -- recall's out of its `pass` string, a sweep's out
/// of the run's own score column. So it works over the frame, which is given
/// the verdict its format calls for, rather than over either reader.
pub fn summary(path: &Path, out: &Path) -> anyhow::Result<()> {
    let mut scores = Frame::open(path)?;

    ensure!(
        scores.format() == runs::FORMAT,
        "{} opens `#= format {}`, which is not a runs table",
        path.display(),
        scores.format()
    );

    scores.layout(runs::layout(&scores.meta));
    scores.verdict(Verdict {
        at: runs::SCORES_AT,
        cutoffs: Named::read(&scores.meta.cutoffs, scores.meta.c)?,
    });

    let hmmer = scores.meta.hmmer()?;
    let runs = scores.meta.runs.len();

    let mut found = vec![0usize; runs];
    let mut hits = vec![0usize; runs];
    let mut hits_sd = vec![0usize; runs];
    let (mut truth, mut truth_sd) = (0usize, 0usize);
    let mut rows = 0u64;

    while scores.step()? {
        rows += 1;

        let true_hit = scores.passed(hmmer);

        // a hit hmmer breaks into one region is a different question from one
        // it breaks into several: the tools disagree most about the second
        let single = true_hit && scores.domain_count() == 1;

        if true_hit {
            truth += 1;
            truth_sd += usize::from(single);
        }

        for run in 0..runs {
            if !scores.passed(run) {
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
    }

    ensure!(
        truth > 0,
        "hmmer found nothing that clears a cutoff; there is nothing to measure against"
    );

    // every setting any run recorded, so a pipeline that swept two knobs gets
    // two columns and one that swept none gets none
    let keys: BTreeSet<&str> = scores
        .meta
        .runs
        .iter()
        .flat_map(|run| run.params.keys())
        .map(String::as_str)
        .collect();

    let mut headers = vec!["name".to_string(), "tool".to_string()];
    headers.extend(keys.iter().map(|k| k.to_string()));
    headers.extend(
        ["wall_s", "found", "hits", "sens", "hits_sd", "sens_sd"]
            .iter()
            .map(|h| h.to_string()),
    );

    let cells: Vec<Vec<String>> = scores
        .meta
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
            ]);

            cells
        })
        .collect();

    let mut table = toil::Table::new(toil::Schema::new(headers));
    preamble(&mut table, &scores.meta, truth, rows);
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
}

impl Tally {
    fn new(runs: usize) -> Tally {
        Tally {
            truth: 0,
            lost_seed: vec![0; runs],
            lost_cloud_align: vec![0; runs],
            lost_cutoff: vec![0; runs],
            kept: vec![0; runs],
        }
    }
}

/// Where the hits hmmer found are lost, per unit and per run.
pub fn stages(path: &Path, out: &Path) -> anyhow::Result<()> {
    let mut scores = Runs::open(path)?;

    let hmmer = scores.meta().hmmer()?;
    let runs = scores.meta().runs.len();

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

    while scores.step()? {
        rows += 1;

        if !scores.row().passed(hmmer) {
            continue;
        }

        let unit = scores.row().shard().to_string();
        let tally = at.entry(unit).or_insert_with(|| Tally::new(runs));
        tally.truth += 1;

        for run in 0..runs {
            if run == hmmer {
                continue;
            }

            // per run rather than per pair: a seeding sweep gives every arm
            // its own seed list, so whether the pair was ever offered is the
            // arm's answer and not the pipeline's
            match (
                scores.seeded(run),
                scores.present(run),
                scores.row().passed(run),
            ) {
                (false, _, _) => tally.lost_seed[run] += 1,
                (_, false, _) => tally.lost_cloud_align[run] += 1,
                (_, _, false) => tally.lost_cutoff[run] += 1,
                (_, _, true) => tally.kept[run] += 1,
            }
        }
    }

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
            ]);
        }
    }

    ensure!(
        !cells.is_empty(),
        "nothing but hmmer ran, so there is no pipeline to trace"
    );

    let truth: usize = at.values().map(|tally| tally.truth).sum();

    let mut table = toil::Table::new(toil::Schema::new(headers));
    preamble(&mut table, scores.meta(), truth, rows);
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
