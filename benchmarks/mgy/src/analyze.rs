//! The analyses, which are all groupings over scores.tbl.
//!
//! Nothing in here reads a results file or knows which pipeline produced the
//! table. A summary is a count per column, a funnel is a count per checkpoint;
//! both are held to the same denominator, which is what hmmer found and scored
//! over its family's cutoff.
//!
//! They are separate from `parse` because reading every results table is the
//! expensive half and the half least likely to change: a different statistic
//! is a re-run of this, not of the benchmark.

use std::collections::BTreeSet;
use std::path::Path;

use anyhow::{Context, ensure};

use crate::scores::{Scores, Tool};

/// What every run found, and what it cost.
///
/// One row per column of scores.tbl, with each run's own settings carried
/// across as columns -- so the same table is cloud-search's (A, B) surface,
/// recall's sensitivity per prefilter setting, and the recall-against-runtime
/// points, depending only on what the pipeline swept.
pub fn summary(scores: &Scores, out: &Path) -> anyhow::Result<()> {
    let hmmer = scores.hmmer()?;

    let truth = scores.denominator(hmmer);
    ensure!(
        truth > 0,
        "hmmer found nothing that clears a cutoff; there is nothing to measure against"
    );

    // a hit hmmer breaks into one region is a different question from one it
    // breaks into several: the tools disagree most about the second kind
    let single: Vec<bool> = scores
        .rows
        .iter()
        .map(|r| r.domain_count(hmmer) == 1)
        .collect();

    let truth_sd = scores
        .rows
        .iter()
        .zip(&single)
        .filter(|(r, sd)| **sd && r.clears(Tool::Hmmer, r.scores[hmmer]))
        .count();

    // every setting any run recorded, so a pipeline that swept two knobs gets
    // two columns and one that swept none gets none
    let keys: BTreeSet<&str> = scores
        .runs
        .iter()
        .flat_map(|r| r.params.keys())
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
        .runs
        .iter()
        .enumerate()
        .map(|(i, run)| {
            let mut found = 0usize;
            let mut hits = 0usize;
            let mut hits_sd = 0usize;

            for (row, &sd) in scores.rows.iter().zip(&single) {
                if !row.clears(run.tool, row.scores[i]) {
                    continue;
                }

                found += 1;

                // held to hmmer as well, so a run is credited for what it
                // agreed with rather than for everything it scored highly
                if row.clears(Tool::Hmmer, row.scores[hmmer]) {
                    hits += 1;
                    hits_sd += usize::from(sd);
                }
            }

            let mut cells = vec![run.name.clone(), run.tool.to_string()];
            cells.extend(keys.iter().map(|k| match run.params.get(*k) {
                Some(value) => value.clone(),
                None => "-".to_string(),
            }));
            cells.extend([
                format!("{:.4}", run.wall_s),
                found.to_string(),
                hits.to_string(),
                format!("{:.4}", frac(hits, truth)),
                hits_sd.to_string(),
                format!("{:.4}", frac(hits_sd, truth_sd)),
            ]);

            cells
        })
        .collect();

    write(out, &meta(scores, hmmer, truth), &headers, &cells)
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

    write(out, &meta(scores, hmmer, truth), &headers, &cells)
}

fn frac(n: usize, of: usize) -> f64 {
    match of {
        0 => 0.0,
        of => n as f64 / of as f64,
    }
}

/// What was searched, what the fractions are fractions of, and the two times
/// the figures want as reference lines.
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

/// A padded table under a `#` header, the way every other table here is
/// written: the rule under the header is dashes, and a row's cells sit under
/// the names rather than beside them.
fn write(path: &Path, meta: &str, headers: &[String], cells: &[Vec<String>]) -> anyhow::Result<()> {
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

    let mut out = meta.to_string();

    out.push('#');
    for (h, &w) in headers.iter().zip(&widths) {
        out.push_str(&format!(" {h:<w$}"));
    }

    out.push_str("\n#");
    for &w in &widths {
        out.push_str(&format!(" {}", "-".repeat(w)));
    }
    out.push('\n');

    for row in cells {
        // the two the `# ` takes on a header line
        out.push_str("  ");
        for (i, (c, &w)) in row.iter().zip(&widths).enumerate() {
            if i > 0 {
                out.push(' ');
            }
            out.push_str(&format!("{c:<w$}"));
        }
        out.push('\n');
    }

    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).context("failed to make the output directory")?;
    }

    std::fs::write(path, out).with_context(|| format!("failed to write {}", path.display()))
}
