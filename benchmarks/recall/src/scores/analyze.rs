//! The analyses, which are groupings over a table of pairs.
//!
//! Nothing in here reads a results file. A summary counts per run, and stages
//! counts per checkpoint within a run; both are held to the same denominator,
//! which is what hmmer found and scored over its family's cutoff.
//!
//! They are separate from `parse` because reading every results table is the
//! expensive half and the half least likely to change: a different statistic
//! is a re-run of this, not of the benchmark.

use std::collections::{BTreeSet, HashSet};
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
pub fn summary(path: &Path, out: &Path, threads: usize) -> anyhow::Result<()> {
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

    let parts = scores.par_each(
        path,
        threads,
        || Tally::new(runs),
        |t, row| {
            t.rows += 1;

            let true_hit = row.passed(hmmer);
            let by_domain = row.domain_passed(hmmer);

            // a hit one of hmmer's domains carries on its own is a different
            // question from one only the whole sequence adds up to: the tools
            // disagree most about the second
            let dom = true_hit && by_domain;
            let seq = true_hit && !by_domain;

            if true_hit {
                t.truth += 1;
                t.truth_dom += usize::from(dom);
                t.truth_seq += usize::from(seq);
            } else if by_domain {
                t.domain_only += 1;
            }

            for run in 0..runs {
                if !row.passed(run) {
                    continue;
                }

                t.found[run] += 1;

                // held to hmmer as well, so a run is credited for what it
                // agreed with rather than for everything it scored highly
                if true_hit {
                    t.hits[run] += 1;
                    t.hits_dom[run] += usize::from(dom);
                    t.hits_seq[run] += usize::from(seq);
                }
            }

            Ok(())
        },
    )?;

    let Tally {
        found,
        hits,
        hits_dom,
        hits_seq,
        truth,
        truth_dom,
        truth_seq,
        domain_only,
        rows,
    } = parts.into_iter().fold(Tally::new(runs), Tally::add);

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

/// What one part of a summary's pass counted.
struct Tally {
    found: Vec<usize>,
    hits: Vec<usize>,
    hits_dom: Vec<usize>,
    hits_seq: Vec<usize>,
    truth: usize,
    truth_dom: usize,
    truth_seq: usize,
    /// Pairs whose best domain cleared while the sequence score did not,
    /// which are no part of the truth and are counted to say how rare.
    domain_only: usize,
    rows: u64,
}

impl Tally {
    fn new(runs: usize) -> Tally {
        Tally {
            found: vec![0; runs],
            hits: vec![0; runs],
            hits_dom: vec![0; runs],
            hits_seq: vec![0; runs],
            truth: 0,
            truth_dom: 0,
            truth_seq: 0,
            domain_only: 0,
            rows: 0,
        }
    }

    fn add(mut self, other: Tally) -> Tally {
        for (mine, theirs) in [
            (&mut self.found, &other.found),
            (&mut self.hits, &other.hits),
            (&mut self.hits_dom, &other.hits_dom),
            (&mut self.hits_seq, &other.hits_seq),
        ] {
            for (a, b) in mine.iter_mut().zip(theirs) {
                *a += b;
            }
        }

        self.truth += other.truth;
        self.truth_dom += other.truth_dom;
        self.truth_seq += other.truth_seq;
        self.domain_only += other.domain_only;
        self.rows += other.rows;
        self
    }
}

/// What every run reported over its own cutoff, split by what hmmer made of
/// the same pairs: its best domain over `hmmer_dom`, reported with no domain
/// over it, or not reported at all.
///
/// Nothing here is a fraction of hmmer's hits. hmmer counts a pair when its
/// sequence score clears `hmmer` or its best domain clears `hmmer_dom`. A
/// family short of any cutoff is left out of every row, and `#= left_out`
/// counts it and the pairs each tool reported in it.
pub fn coverage(path: &Path, out: &Path, threads: usize) -> anyhow::Result<()> {
    let mut scores = Frame::open(path)?;

    ensure!(
        scores.format() == crate::scores::FORMAT,
        "{} opens `#= format {}`, which is not recall's table",
        path.display(),
        scores.format()
    );

    scores.layout(read::layout(&scores.meta));

    let hmmer = scores.meta.hmmer()?;
    let hmmer_score = read::score(&scores.meta, Tool::Hmmer).context("no hmmer score column")?;
    let complete = crate::scores::complete(&scores.meta.cutoffs, scores.meta.c)?;

    let tools = scores.meta.tools();
    let at: Vec<usize> = tools
        .iter()
        .map(|tool| read::score(&scores.meta, *tool).context("a tool with no score column"))
        .collect::<anyhow::Result<_>>()?;

    let runs = scores.meta.runs.len();

    let parts = scores.par_each(
        path,
        threads,
        || Cover {
            counts: vec![[0; 3]; runs],
            left_families: HashSet::new(),
            left: vec![0; at.len()],
        },
        |t, row| {
            let query =
                std::str::from_utf8(row.field(0)).context("a query name that is not utf-8")?;

            if !complete.contains(query) {
                if !t.left_families.contains(query) {
                    t.left_families.insert(query.to_string());
                }
                for (i, at) in at.iter().enumerate() {
                    t.left[i] += usize::from(row.field(*at) != b"-");
                }
                return Ok(());
            }

            let by_domain = row.domain_passed(hmmer);
            let class = match (by_domain, row.field(hmmer_score) != b"-") {
                (true, _) => 0,
                (false, true) => 1,
                (false, false) => 2,
            };

            for run in 0..runs {
                if row.passed(run) || (run == hmmer && by_domain) {
                    t.counts[run][class] += 1;
                }
            }

            Ok(())
        },
    )?;

    let mut counts = vec![[0usize; 3]; runs];
    let mut left_families = HashSet::new();
    let mut left = vec![0usize; tools.len()];
    for part in parts {
        for (mine, theirs) in counts.iter_mut().zip(&part.counts) {
            for (a, b) in mine.iter_mut().zip(theirs) {
                *a += b;
            }
        }
        for (a, b) in left.iter_mut().zip(&part.left) {
            *a += b;
        }
        left_families.extend(part.left_families);
    }

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
        ["wall_s", "hits", "dom_over", "dom_under", "no_hmmer"]
            .iter()
            .map(|h| h.to_string()),
    );

    let mut table = toil::Table::new(toil::Schema::new(headers));
    table.meta(
        "left_out",
        [toil::Cell::from("families"), left_families.len().into()],
    );
    table.meta(
        "left_out",
        tools
            .iter()
            .zip(&left)
            .flat_map(|(tool, n)| [toil::Cell::from(tool.to_string()), (*n).into()]),
    );

    for (run, count) in scores.meta.runs.iter().zip(&counts) {
        let mut cells = vec![run.name.clone(), run.tool.to_string()];
        cells.extend(keys.iter().map(|k| match run.params.get(*k) {
            Some(value) => value.clone(),
            None => "-".to_string(),
        }));
        cells.push(format!("{:.4}", run.wall_s));
        cells.push(count.iter().sum::<usize>().to_string());
        cells.extend(count.iter().map(usize::to_string));

        table.row(cells);
    }

    table
        .write(out)
        .with_context(|| format!("failed to write {}", out.display()))
}

/// What one part of a coverage pass counted.
struct Cover {
    /// Per run: dom_over, dom_under, no_hmmer.
    counts: Vec<[usize; 3]>,
    left_families: HashSet<String>,
    /// Per tool, the pairs it reported in a family left out.
    left: Vec<usize>,
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
