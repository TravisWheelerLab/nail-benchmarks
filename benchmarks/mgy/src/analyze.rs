//! The analyses, which are groupings over a table of pairs.
//!
//! Nothing in here reads a results file. A summary is a count per run, a
//! funnel is a count per checkpoint; both are held to the same denominator,
//! which is what hmmer found and scored over its family's cutoff.
//!
//! They are separate from `parse` because reading every results table is the
//! expensive half and the half least likely to change: a different statistic
//! is a re-run of this, not of the benchmark.

use std::collections::BTreeSet;
use std::path::Path;

use anyhow::ensure;

use util::tbl;

use crate::scores::read::Reader;
use crate::scores::{self, Scores, Tool};

/// What every run found, and what it cost.
///
/// One pass over `scores.tbl`, counting per run: what it reported over the
/// cutoff, how much of that hmmer also found, and how much of that hmmer read
/// as one domain rather than several. Nothing is held but the counters.
pub fn summary(path: &Path, out: &Path) -> anyhow::Result<()> {
    let mut scores = Reader::open(path)?;

    let hmmer = scores.meta.hmmer()?;
    let runs = scores.meta.runs.len();

    let mut found = vec![0usize; runs];
    let mut hits = vec![0usize; runs];
    let mut hits_sd = vec![0usize; runs];
    let (mut truth, mut truth_sd) = (0usize, 0usize);
    let mut rows = 0u64;

    while let Some(row) = scores.next()? {
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

    tbl::write(
        out,
        tbl::Table {
            meta: &preamble(&scores.meta, truth, rows),
            headers: &headers,
            rows: &cells,
            ragged_last: false,
        },
    )
}

/// What was searched, what the fractions are fractions of, and the two times
/// the figures use as reference lines.
fn preamble(meta: &scores::Meta, truth: usize, rows: u64) -> String {
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

    format!(
        "# query  {:>9} families  {:>12} residues  {:>12} bytes\n\
         # target {:>9} seqs      {:>12} residues  {:>12} bytes\n\
         # pairs  {:>9} rows      {:>12} runs\n\
         # hmmer  {:>9} hits      {:>12.4} wall_s\n\
         # seed   {:>9}           {:>12} wall_s\n\
         #\n",
        meta.query.count,
        meta.query.residues,
        meta.query.bytes,
        count,
        residues,
        bytes,
        rows,
        meta.runs.len(),
        truth,
        hmmer,
        "",
        // recall does not seed: nail's prefilter is part of its search, and
        // there is no stage of its own to time
        "-",
    )
}

/// Where the hits hmmer found are lost, for every run that isn't hmmer's.
///
/// Only two checkpoints are visible from the outside: whether a pair got a
/// seed, and whether it ended up in the tool's table at all. Everything
/// between them collapses into one bucket -- see hit_loss.rs for why the
/// e-value gate has to be opened up for that bucket to mean what it says.
///
/// Reaching the table is presence, not a cutoff: this is asking where a pair
/// was dropped, and a pair that survived to be scored badly was not dropped.
pub fn funnel(scores: &Scores, out: &Path) -> anyhow::Result<()> {
    let hmmer = scores.hmmer()?;

    ensure!(
        scores.rows.iter().any(|r| r.seeded.is_some()),
        "this pipeline kept no seeds, so there is no seeding checkpoint to split on"
    );

    let truth = scores.denominator(hmmer);
    ensure!(
        truth > 0,
        "hmmer found nothing that clears a cutoff; there is nothing to measure against"
    );

    let headers = ["run", "stage", "n", "sens"].map(str::to_string).to_vec();
    let mut cells: Vec<Vec<String>> = Vec::new();

    for (i, run) in scores.runs.iter().enumerate() {
        if i == hmmer {
            continue;
        }

        let (mut lost_seed, mut lost_cloud_align, mut reported) = (0usize, 0usize, 0usize);

        for row in &scores.rows {
            if !row.clears(Tool::Hmmer, row.scores[hmmer]) {
                continue;
            }

            match (row.seeded, row.scores[i]) {
                (Some(false), _) => lost_seed += 1,
                (_, None) => lost_cloud_align += 1,
                (_, Some(_)) => reported += 1,
            }
        }

        // each stage is what it dropped, against what is still standing after
        // it -- so the last column falls from 1 to the fraction that survived
        let stages = [
            ("truth", truth, truth),
            ("lost_seed", lost_seed, truth - lost_seed),
            (
                "lost_cloud_align",
                lost_cloud_align,
                truth - lost_seed - lost_cloud_align,
            ),
            ("reported", reported, reported),
        ];

        cells.extend(stages.iter().map(|&(stage, n, left)| {
            vec![
                run.name.clone(),
                stage.to_string(),
                n.to_string(),
                format!("{:.4}", frac(left, truth)),
            ]
        }));
    }

    ensure!(
        !cells.is_empty(),
        "nothing but hmmer ran, so there is no pipeline to trace"
    );

    tbl::write(
        out,
        tbl::Table {
            meta: &meta(scores, hmmer, truth),
            headers: &headers,
            rows: &cells,
            ragged_last: false,
        },
    )
}

fn frac(n: usize, of: usize) -> f64 {
    match of {
        0 => 0.0,
        of => n as f64 / of as f64,
    }
}

/// What was searched, what the fractions are fractions of, and the two times
/// the figures use as reference lines.
fn meta(scores: &Scores, hmmer: usize, truth: usize) -> String {
    let (mut count, mut residues, mut bytes) = (0usize, 0u64, 0u64);
    for (_, size) in &scores.targets {
        count += size.count;
        residues += size.residues;
        bytes += size.bytes;
    }

    // a dash rather than a zero for the pipelines that never seeded: seeding
    // taking no time and there being no seeding are different things
    let seed_wall_s = match scores.seed_wall_s.is_empty() {
        true => "-".to_string(),
        false => format!(
            "{:.4}",
            scores.seed_wall_s.iter().map(|(_, w)| w).sum::<f64>()
        ),
    };

    format!(
        "# query  {:>9} families  {:>12} residues  {:>12} bytes\n\
         # target {:>9} seqs      {:>12} residues  {:>12} bytes\n\
         # pairs  {:>9} rows      {:>12} runs\n\
         # hmmer  {:>9} hits      {:>12.4} wall_s\n\
         # seed   {:>9}           {:>12} wall_s\n\
         #\n",
        scores.query.count,
        scores.query.residues,
        scores.query.bytes,
        count,
        residues,
        bytes,
        scores.rows.len(),
        scores.runs.len(),
        truth,
        scores.runs[hmmer].wall_s,
        "",
        seed_wall_s,
    )
}
