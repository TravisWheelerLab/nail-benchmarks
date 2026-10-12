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

use crate::scores::frame::Frame;
use crate::scores::{Tool, read};

/// What every run found, and what it cost.
///
/// One pass over the table, counting per run: what it reported over the
/// cutoff, how much of that hmmer also found, and that split by whether one
/// of hmmer's domains cleared its own cutoff or only the sequence score did.
/// Nothing is held but the counters.
///
/// This reads either table. A summary asks only whether each run cleared its
/// family's cutoff and what hmmer's domain letter says, and the two tables answer
/// the first differently -- recall's out of its `pass` string, a sweep's out
/// of the run's own score column. So it works over the frame, which is given
/// the verdict its format calls for, rather than over either reader.
pub fn summary(path: &Path, out: &Path) -> anyhow::Result<()> {
    let mut scores = Frame::open(path)?;

    ensure!(
        scores.format() == crate::scores::FORMAT,
        "{} opens `#= format {}`, which is not recall's table",
        path.display(),
        scores.format()
    );

    scores.layout(read::layout(&scores.meta));

    let hmmer = scores.meta.hmmer()?;
    let runs = scores.meta.runs.len();

    let mut found = vec![0usize; runs];
    let mut hits = vec![0usize; runs];
    let mut hits_dom = vec![0usize; runs];
    let mut hits_seq = vec![0usize; runs];
    let (mut truth, mut truth_dom, mut truth_seq) = (0usize, 0usize, 0usize);
    // pairs whose best domain cleared while the sequence score did not,
    // which are no part of the truth and are counted to say how rare
    let mut domain_only = 0usize;
    let mut rows = 0u64;

    scores.each(|row| {
        rows += 1;

        let true_hit = row.passed(hmmer);
        let by_domain = row.domain_passed(hmmer);

        // a hit one of hmmer's domains carries on its own is a different
        // question from one only the whole sequence adds up to: the tools
        // disagree most about the second
        let dom = true_hit && by_domain;
        let seq = true_hit && !by_domain;

        if true_hit {
            truth += 1;
            truth_dom += usize::from(dom);
            truth_seq += usize::from(seq);
        } else if by_domain {
            domain_only += 1;
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
                hits_dom[run] += usize::from(dom);
                hits_seq[run] += usize::from(seq);
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
        [
            "wall_s", "found", "hits", "sens", "hits_dom", "sens_dom", "hits_seq", "sens_seq",
        ]
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
                hits_dom[i].to_string(),
                format!("{:.4}", frac(hits_dom[i], truth_dom)),
                hits_seq[i].to_string(),
                format!("{:.4}", frac(hits_seq[i], truth_seq)),
            ]);

            cells
        })
        .collect();

    let mut table = toil::Table::new(toil::Schema::new(headers));
    preamble(
        &mut table,
        &scores.meta,
        [truth, truth_dom, truth_seq, domain_only],
        rows,
    );
    for row in cells {
        table.row(row);
    }

    table
        .write(out)
        .with_context(|| format!("failed to write {}", out.display()))
}

/// What was searched, what the fractions are fractions of, and the two times
/// the figures use as reference lines.
///
/// `truth` is hmmer's hits, then how many a domain carried, how many only the
/// sequence score did, and the pairs where a domain cleared and the sequence
/// score did not, which are in none of the fractions.
fn preamble(table: &mut toil::Table, meta: &crate::scores::Meta, truth: [usize; 4], rows: u64) {
    let [truth, dom, seq, domain_only] = truth;
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
        .comment(format!(
            "       {dom:>9} by domain {seq:>12} by sequence only {domain_only:>12} domain only"
        ))
        .comment(format!("seed   {:>9}           {seed:>12} wall_s", ""))
        .comment("");
}

fn frac(n: usize, of: usize) -> f64 {
    match of {
        0 => 0.0,
        of => n as f64 / of as f64,
    }
}
