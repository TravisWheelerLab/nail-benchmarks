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

use anyhow::{bail, ensure};

use util::tbl;

use crate::frame::Frame;
use crate::runs::Reader as Runs;
use crate::{Tool, read, runs};

/// What every run found, and what it cost.
///
/// One pass over the table, counting per run: what it reported over the
/// cutoff, how much of that hmmer also found, and how much of that hmmer read
/// as one domain rather than several. Nothing is held but the counters.
///
/// This reads either table. A summary asks only what the `pass` string and
/// the domain list say, and both tables carry those in the same place -- what
/// tells them apart is the score columns, which a summary never opens. So it
/// works over the frame rather than over either reader.
pub fn summary(path: &Path, out: &Path) -> anyhow::Result<()> {
    let mut scores = Frame::open(path)?;

    let layout = match scores.format() {
        f if f == crate::FORMAT => read::layout(&scores.meta),
        f if f == runs::FORMAT => runs::layout(&scores.meta),
        other => bail!("{} opens `{other}`, which is no table here", path.display()),
    };
    scores.layout(layout);

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
fn preamble(meta: &crate::Meta, truth: usize, rows: u64) -> String {
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
        seed,
    )
}

/// Where the hits hmmer found are lost, for every run that isn't hmmer's.
///
/// Only two checkpoints are visible from the outside: whether a pair got a
/// seed, and whether it ended up in the run's table at all. Everything
/// between them collapses into one bucket -- see loss-decomp for why the
/// e-value gate has to be opened up for that bucket to mean what it says.
///
/// Reaching the table is presence, not a cutoff: this is asking where a pair
/// was dropped, and a pair that survived to be scored badly was not dropped.
/// That is why this reads `runs.tbl` rather than recall's table -- presence is
/// a property of a run, and recall keeps a column per tool.
///
/// One pass, counting per run. A sweep that seeded once and searched every
/// cell off those seeds gets a column per cell, which is what tells the loss
/// every cell shares from the loss its pruning caused.
/// What one unit's pairs came to, per run.
struct Tally {
    /// Pairs hmmer reported at or above their family's cutoff, for this unit.
    truth: usize,
    lost_seed: Vec<usize>,
    lost_cloud_align: Vec<usize>,
    reported: Vec<usize>,
}

impl Tally {
    fn new(runs: usize) -> Tally {
        Tally {
            truth: 0,
            lost_seed: vec![0; runs],
            lost_cloud_align: vec![0; runs],
            reported: vec![0; runs],
        }
    }
}

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
            match (scores.seeded(run), scores.present(run)) {
                (false, _) => tally.lost_seed[run] += 1,
                (_, false) => tally.lost_cloud_align[run] += 1,
                (_, true) => tally.reported[run] += 1,
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
    let headers = ["unit", "run", "truth", "lost_seed", "lost_align", "reported", "sens"]
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
                tally.reported[run].to_string(),
                format!("{:.4}", frac(tally.reported[run], tally.truth)),
            ]);
        }
    }

    ensure!(
        !cells.is_empty(),
        "nothing but hmmer ran, so there is no pipeline to trace"
    );

    let truth: usize = at.values().map(|tally| tally.truth).sum();

    tbl::write(
        out,
        tbl::Table {
            meta: &preamble(scores.meta(), truth, rows),
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
