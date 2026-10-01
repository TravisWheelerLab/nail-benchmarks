//! The analyses, which are groupings over a table of pairs.
//!
//! `stages` counts per checkpoint within a run, held to one denominator:
//! what hmmer found and scored over its family's cutoff. It reads no results
//! table, and reads each arm's prefilter database for one thing: whether a
//! pair the seeding never offered was in the prefilter list at all.
//!
//! They are separate from `parse` because reading every results table is the
//! expensive half and the half least likely to change: a different statistic
//! is a re-run of this, not of the benchmark.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use anyhow::{Context, ensure};

use crate::scores::Tool;
use crate::scores::depth;
use crate::scores::runs::Reader as Runs;

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
    /// Not in the arm's prefilter list at any rank.
    lost_prefilter: Vec<usize>,
    /// In the prefilter list, and not in the seed list: mmseqs' alignment
    /// dropped it.
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
            lost_prefilter: vec![0; runs],
            lost_seed: vec![0; runs],
            lost_cloud_align: vec![0; runs],
            lost_cutoff: vec![0; runs],
            kept: vec![0; runs],
        }
    }
}

/// Which of hmmer's hits a table is held to, by how many domains hmmer
/// resolved each into.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Domains {
    All,
    Single,
    Multi,
}

impl Domains {
    const EACH: [Domains; 3] = [Domains::All, Domains::Single, Domains::Multi];

    fn holds(self, single: bool) -> bool {
        match self {
            Domains::All => true,
            Domains::Single => single,
            Domains::Multi => !single,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Domains::All => "all",
            Domains::Single => "single",
            Domains::Multi => "multi",
        }
    }

    /// Where this set's table goes, beside the one `out` names: `x.tbl`,
    /// `x-single.tbl`, `x-multi.tbl`.
    fn path(self, out: &Path) -> PathBuf {
        match self {
            Domains::All => out.to_path_buf(),
            set => {
                let stem = out.file_stem().unwrap_or_default().to_string_lossy();
                let ext = out
                    .extension()
                    .map(|e| e.to_string_lossy())
                    .unwrap_or_default();
                out.with_file_name(format!("{stem}-{}.{ext}", set.name()))
            }
        }
    }
}

/// Where the hits hmmer found are lost, per unit and per run, three times
/// over: for every hit, for the hits hmmer resolved as one domain, and for
/// the rest.
///
/// `root` is the run directory, for the prefilter databases under its
/// `results/`. `out` names the table of every hit, and the other two go
/// beside it.
pub fn stages(path: &Path, root: &Path, out: &Path) -> anyhow::Result<()> {
    let mut scores = Runs::open(path)?;

    let hmmer = scores.meta().hmmer()?;
    let runs = scores.meta().runs.len();

    ensure!(
        scores.meta().runs.iter().any(|run| run.seeds.is_some()),
        "no run kept its seed list, so there is no seeding checkpoint to split on"
    );

    // per unit as well as per run. A `cross` set searches one query against
    // several kinds of target, and what hmmer found in one is not the truth
    // set for another: summed, the two corpora make a sensitivity that
    // describes neither
    let mut at: indexmap::IndexMap<String, [Tally; 3]> = indexmap::IndexMap::new();
    let mut rows = 0u64;

    // the pairs no seed list held, per unit and run, with whether hmmer
    // resolved each as one domain, to be looked up in the arm's prefilter
    // list once the table has been read
    let mut unseeded: HashMap<(String, usize), Vec<(String, String, bool)>> = HashMap::new();

    scores.each(|row| {
        rows += 1;

        if !row.row().passed(hmmer) {
            return Ok(());
        }

        let unit = row.row().shard().to_string();
        let single = row.domain_count() == 1;
        let tallies = at
            .entry(unit.clone())
            .or_insert_with(|| std::array::from_fn(|_| Tally::new(runs)));

        for (set, tally) in Domains::EACH.iter().zip(tallies.iter_mut()) {
            if !set.holds(single) {
                continue;
            }
            tally.truth += 1;

            for run in 0..runs {
                if run == hmmer {
                    continue;
                }

                // per run rather than per pair: a seeding sweep gives every
                // arm its own seed list, so whether the pair was ever offered
                // is the arm's answer and not the pipeline's
                match (row.seeded(run), row.present(run), row.row().passed(run)) {
                    (false, _, _) => tally.lost_seed[run] += 1,
                    (_, false, _) => tally.lost_cloud_align[run] += 1,
                    (_, _, false) => tally.lost_cutoff[run] += 1,
                    (_, _, true) => tally.kept[run] += 1,
                }
            }
        }

        for run in 0..runs {
            if run != hmmer && !row.seeded(run) {
                unseeded.entry((unit.clone(), run)).or_default().push((
                    String::from_utf8_lossy(row.row().field(0)).into_owned(),
                    String::from_utf8_lossy(row.row().field(1)).into_owned(),
                    single,
                ));
            }
        }

        Ok(())
    })?;

    ensure!(
        at.values().any(|t| t[0].truth > 0),
        "hmmer found nothing that clears a cutoff; there is nothing to measure against"
    );

    // a pair the prefilter never returned could not have been seeded, so
    // the pairs no seed list held split on the prefilter list: not there at
    // all, or there and dropped by mmseqs' alignment
    for ((unit, run), pairs) in unseeded {
        let column = &scores.meta().runs[run];
        let arm = column
            .seeds
            .as_deref()
            .with_context(|| format!("run {:?} kept no seed list", column.name))?;

        let mut names: Vec<(String, String)> = pairs
            .iter()
            .map(|(q, t, _)| (q.clone(), t.clone()))
            .collect();
        let missing: HashSet<(String, String)> = depth::beyond(root, arm, &unit, &mut names)?
            .into_iter()
            .collect();

        let tallies = at
            .get_mut(&unit)
            .expect("every unit with pairs was tallied");
        for (q, t, single) in pairs {
            if !missing.contains(&(q, t)) {
                continue;
            }
            for (set, tally) in Domains::EACH.iter().zip(tallies.iter_mut()) {
                if set.holds(single) {
                    tally.lost_prefilter[run] += 1;
                    tally.lost_seed[run] -= 1;
                }
            }
        }
    }

    // one row per (unit, run), a column per checkpoint. Written long it was
    // four rows apiece, where `n` meant a population on two of them and a loss
    // on the other two, and the last row's fraction only ever repeated the one
    // above it
    let headers = [
        "unit",
        "run",
        "truth",
        "lost_prefilter",
        "lost_seed",
        "lost_align",
        "lost_cutoff",
        "kept",
        "sens",
    ]
    .map(str::to_string)
    .to_vec();

    for (i, set) in Domains::EACH.iter().enumerate() {
        let mut cells: Vec<Vec<String>> = Vec::new();

        for (unit, tallies) in &at {
            let tally = &tallies[i];
            for (run, column) in scores.meta().runs.iter().enumerate() {
                if run == hmmer {
                    continue;
                }

                cells.push(vec![
                    unit.clone(),
                    column.name.clone(),
                    tally.truth.to_string(),
                    tally.lost_prefilter[run].to_string(),
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

        let truth: usize = at.values().map(|tallies| tallies[i].truth).sum();

        let mut table = toil::Table::new(toil::Schema::new(headers.clone()));
        preamble(&mut table, scores.meta(), truth, rows);
        table.meta("domains", [set.name()]);
        for row in cells {
            table.row(row);
        }

        let path = set.path(out);
        table
            .write(&path)
            .with_context(|| format!("failed to write {}", path.display()))?;
    }

    Ok(())
}

fn frac(n: usize, of: usize) -> f64 {
    match of {
        0 => 0.0,
        of => n as f64 / of as f64,
    }
}
