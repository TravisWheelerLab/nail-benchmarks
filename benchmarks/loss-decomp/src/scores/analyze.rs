//! The analyses, which are groupings over a table of pairs.
//!
//! `stages` counts per checkpoint within a run, held to one denominator:
//! what hmmer found and scored over its family's cutoff. It reads no results
//! table, and reads each arm's prefilter database for one thing: whether a
//! pair the seeding never offered was in the prefilter list at all.
//!
//! A family short of a cutoff for some tool is left out of every count, and
//! `#= left_out` counts the families and the pairs each run reported in them.
//!
//! They are separate from `parse` because reading every results table is the
//! expensive half and the half least likely to change: a different statistic
//! is a re-run of this, not of the benchmark.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::{Context, ensure};

use util::domains::{self, Domain};
use util::manifest;

use crate::scores::depth;
use crate::scores::runs::Reader as Runs;
use crate::scores::{LeftOut, Tool};

/// What was searched, what the fractions are fractions of, and the two times
/// the figures use as reference lines.
///
/// `truth` is every pair hmmer's sequence score or best domain carried over
/// its cutoff, and `domain_only` is the part of it only the best domain did.
fn preamble(
    table: &mut toil::Table,
    meta: &crate::scores::Meta,
    truth: usize,
    domain_only: usize,
    rows: u64,
) {
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
        .comment(format!("       {domain_only:>9} domain only"))
        .comment(format!("seed   {:>9}           {seed:>12} wall_s", ""))
        .comment("");
}

/// What one unit's pairs came to, per run.
struct Tally {
    /// Pairs hmmer reported at or above their family's cutoff, for this unit.
    truth: usize,
    /// Not in the arm's prefilter list at any rank.
    lost_prefilter: Vec<usize>,
    /// In the prefilter list, past the rank where prog seeding stopped for
    /// the query, so mmseqs never aligned it.
    lost_stopped: Vec<usize>,
    /// In the prefilter list and aligned by mmseqs, whose alignment dropped
    /// it from the seed list.
    lost_mmseqs: Vec<usize>,
    lost_cloud_align: Vec<usize>,
    /// Scored, and below the family's cutoff. nail found the pair and would
    /// not return it, which is a loss like any other.
    lost_cutoff: Vec<usize>,
    /// Of `lost_cutoff`, the pairs hmmer's best domain clears the domain
    /// cutoff on by itself.
    cutoff_domain: Vec<usize>,
    /// Of `lost_cutoff`, the pairs no domain of hmmer's clears the domain
    /// cutoff: hmmer's sequence score reaches its own, no one domain does.
    cutoff_nodomain: Vec<usize>,
    /// Of `cutoff_nodomain`, the pairs hmmer reported one domain for.
    nodomain_single: Vec<usize>,
    /// Of `cutoff_nodomain`, the pairs hmmer reported more than one domain
    /// for.
    nodomain_multi: Vec<usize>,
    kept: Vec<usize>,
    /// Pairs the run reported at or above the cutoff, whether or not hmmer
    /// found them.
    over: Vec<usize>,
    /// Of `truth`, the pairs only hmmer's best domain carried: it cleared
    /// `hmmer_dom` and the sequence score did not clear `hmmer`.
    domain_only: usize,
    /// Pairs the run reported, at any score, in a family short of a cutoff.
    left_out: Vec<usize>,
}

impl Tally {
    fn new(runs: usize) -> Tally {
        Tally {
            truth: 0,
            lost_prefilter: vec![0; runs],
            lost_stopped: vec![0; runs],
            lost_mmseqs: vec![0; runs],
            lost_cloud_align: vec![0; runs],
            lost_cutoff: vec![0; runs],
            cutoff_domain: vec![0; runs],
            cutoff_nodomain: vec![0; runs],
            nodomain_single: vec![0; runs],
            nodomain_multi: vec![0; runs],
            kept: vec![0; runs],
            over: vec![0; runs],
            domain_only: 0,
            left_out: vec![0; runs],
        }
    }
}

/// Where the hits hmmer found are lost, per unit and per run.
///
/// `root` is the run directory, for the prefilter databases under its
/// `results/`. `out` names the table, and `cutoff-pairs.tbl` goes beside it.
pub fn stages(path: &Path, root: &Path, out: &Path, threads: usize) -> anyhow::Result<()> {
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
    let mut at: indexmap::IndexMap<String, Tally> = indexmap::IndexMap::new();
    let mut rows = 0u64;

    // the pairs no seed list held, per unit and run, to be looked up in the
    // arm's prefilter list once the table has been read
    let mut unseeded: HashMap<(String, usize), Vec<(String, String)>> = HashMap::new();

    // every pair a run scored under the cutoff, for the per-pair table:
    // what hmmer and the run scored it, and how hmmer's domains add up
    let mut under: Vec<Under> = Vec::new();

    let mut left = LeftOut::new(runs);

    scores.each(|row| {
        let unit = row.row().shard().to_string();

        if !row.row().complete() {
            let tally = at
                .entry(unit.clone())
                .or_insert_with(|| Tally::new(runs));
            for run in 0..runs {
                tally.left_out[run] += usize::from(row.present(run));
            }
            let family = std::str::from_utf8(row.row().field(0))
                .context("a query name that is not utf-8")?;
            left.pair(family, |run| row.present(run));
            return Ok(());
        }

        rows += 1;

        let best_domain = row.best_domain();
        let dom_cutoff = row.row().domain_cutoff();
        let domain_clears = row.domain_clears();
        let reported = row.domains_reported();
        let tally = at
            .entry(unit.clone())
            .or_insert_with(|| Tally::new(runs));

        for run in 0..runs {
            if row.row().passed(run) {
                tally.over[run] += 1;
            }
        }

        if !row.hmmer_found(hmmer) {
            return Ok(());
        }

        tally.truth += 1;
        tally.domain_only += usize::from(!row.row().passed(hmmer));

        for run in 0..runs {
            if run == hmmer {
                continue;
            }

            // per run rather than per pair: a seeding sweep gives every
            // arm its own seed list, so whether the pair was ever offered
            // is the arm's answer and not the pipeline's
            match (row.seeded(run), row.present(run), row.row().passed(run)) {
                (false, _, _) => tally.lost_mmseqs[run] += 1,
                (_, false, _) => tally.lost_cloud_align[run] += 1,
                (_, _, false) => {
                    tally.lost_cutoff[run] += 1;
                    match (domain_clears, reported) {
                        (true, _) => tally.cutoff_domain[run] += 1,
                        (false, 1) => {
                            tally.cutoff_nodomain[run] += 1;
                            tally.nodomain_single[run] += 1;
                        }
                        (false, _) => {
                            tally.cutoff_nodomain[run] += 1;
                            tally.nodomain_multi[run] += 1;
                        }
                    }
                    if let (Some(cutoff), Some(hmmer_cutoff), Some(best), Some(score), Some(seq)) = (
                        row.row().cutoff(run),
                        row.row().cutoff(hmmer),
                        best_domain,
                        row.score(run),
                        row.score(hmmer),
                    ) {
                        under.push(Under {
                            unit: unit.clone(),
                            run,
                            query: String::from_utf8_lossy(row.row().field(0)).into_owned(),
                            target: String::from_utf8_lossy(row.row().field(1)).into_owned(),
                            cutoff,
                            hmmer_cutoff,
                            dom_cutoff,
                            hmmer: seq,
                            domain: best,
                            nail: score,
                            reported,
                            kind: match (domain_clears, reported) {
                                (true, _) => "domain",
                                (false, 1) => "single",
                                (false, _) => "multi",
                            },
                        });
                    }
                }
                (_, _, true) => tally.kept[run] += 1,
            }
        }

        for run in 0..runs {
            if run != hmmer && !row.seeded(run) {
                unseeded.entry((unit.clone(), run)).or_default().push((
                    String::from_utf8_lossy(row.row().field(0)).into_owned(),
                    String::from_utf8_lossy(row.row().field(1)).into_owned(),
                ));
            }
        }

        Ok(())
    })?;

    ensure!(
        at.values().any(|t| t.truth > 0),
        "hmmer found nothing that clears a cutoff; there is nothing to measure against"
    );

    // the pairs no seed list held split three ways on the prefilter list:
    // not there at all, there past where prog seeding stopped so never
    // aligned, or there, aligned, and dropped by mmseqs' alignment
    for ((unit, run), pairs) in unseeded {
        let column = &scores.meta().runs[run];
        let arm = column
            .seeds
            .as_deref()
            .with_context(|| format!("run {:?} kept no seed list", column.name))?;

        let mut names = pairs.clone();
        let (missing, stopped) = depth::unseeded(root, arm, &unit, &mut names, threads)?;
        let missing: HashSet<(String, String)> = missing.into_iter().collect();
        let stopped: HashSet<(String, String)> = stopped.into_iter().collect();

        let tally = at
            .get_mut(&unit)
            .expect("every unit with pairs was tallied");
        for pair in pairs {
            if missing.contains(&pair) {
                tally.lost_prefilter[run] += 1;
                tally.lost_mmseqs[run] -= 1;
            } else if stopped.contains(&pair) {
                tally.lost_stopped[run] += 1;
                tally.lost_mmseqs[run] -= 1;
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
        "lost_stopped",
        "lost_mmseqs",
        "lost_align",
        "lost_cutoff",
        "cutoff_domain",
        "cutoff_nodomain",
        "nodomain_single",
        "nodomain_multi",
        "kept",
        "over",
        "sens",
        "left_out",
    ]
    .map(str::to_string)
    .to_vec();

    {
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
                    tally.lost_prefilter[run].to_string(),
                    tally.lost_stopped[run].to_string(),
                    tally.lost_mmseqs[run].to_string(),
                    tally.lost_cloud_align[run].to_string(),
                    tally.lost_cutoff[run].to_string(),
                    tally.cutoff_domain[run].to_string(),
                    tally.cutoff_nodomain[run].to_string(),
                    tally.nodomain_single[run].to_string(),
                    tally.nodomain_multi[run].to_string(),
                    tally.kept[run].to_string(),
                    tally.over[run].to_string(),
                    format!("{:.4}", frac(tally.kept[run], tally.truth)),
                    tally.left_out[run].to_string(),
                ]);
            }
        }

        ensure!(
            !cells.is_empty(),
            "nothing but hmmer ran, so there is no pipeline to trace"
        );

        let truth: usize = at.values().map(|tally| tally.truth).sum();
        let domain_only: usize = at.values().map(|tally| tally.domain_only).sum();

        let mut table = toil::Table::new(toil::Schema::new(headers.clone()));
        preamble(&mut table, scores.meta(), truth, domain_only, rows);
        for words in left.lines(&scores.meta().runs) {
            table.meta("left_out", words);
        }
        for row in cells {
            table.row(row);
        }

        table
            .write(out)
            .with_context(|| format!("failed to write {}", out.display()))?;
    }

    // one row per pair a run scored under the cutoff, beside stages.tbl:
    // the scores and the cutoffs each was held to, how many domains hmmer
    // reported, and the paper's class for how they add up to its score,
    // read off hmmer's domtbl for the unit
    let classes = classify(root, &scores.meta().runs[hmmer].name, &under)?;
    let mut table = toil::Table::new(toil::Schema::new(
        [
            "unit",
            "run",
            "query",
            "target",
            "cutoff",
            "hmmer_cutoff",
            "dom_cutoff",
            "hmmer",
            "domain",
            "nail",
            "reported",
            "kind",
            "class",
        ]
        .map(str::to_string)
        .to_vec(),
    ));
    for words in left.lines(&scores.meta().runs) {
        table.meta("left_out", words);
    }
    for (pair, class) in under.iter().zip(classes) {
        table.row(vec![
            pair.unit.clone(),
            scores.meta().runs[pair.run].name.clone(),
            pair.query.clone(),
            pair.target.clone(),
            format!("{:.1}", pair.cutoff),
            format!("{:.1}", pair.hmmer_cutoff),
            match pair.dom_cutoff {
                Some(cutoff) => format!("{cutoff:.1}"),
                None => "-".to_string(),
            },
            format!("{:.1}", pair.hmmer),
            format!("{:.1}", pair.domain),
            format!("{:.1}", pair.nail),
            pair.reported.to_string(),
            pair.kind.to_string(),
            class,
        ]);
    }
    let path = out.with_file_name("cutoff-pairs.tbl");
    table
        .write(&path)
        .with_context(|| format!("failed to write {}", path.display()))?;

    Ok(())
}

/// One pair a run scored under the cutoff.
struct Under {
    unit: String,
    run: usize,
    query: String,
    target: String,
    /// What the run's score was held to.
    cutoff: f32,
    /// What hmmer's sequence score was held to.
    hmmer_cutoff: f32,
    /// What hmmer's best domain was held to, where the family has one.
    dom_cutoff: Option<f32>,
    /// hmmer's sequence score.
    hmmer: f32,
    /// hmmer's best domain score.
    domain: f32,
    /// The run's score.
    nail: f32,
    /// How many domains hmmer reported.
    reported: usize,
    /// `domain` where the best domain clears its cutoff alone, else `single`
    /// or `multi` by how many hmmer reported.
    kind: &'static str,
}

/// The paper's class for each pair's domains against its cutoff, from the
/// domtbl hmmer wrote for the pair's unit: `strong_match`, a `pileup`,
/// `colinear_match`, `shuffled_match` or `weak_match`, and `-` for a pair the
/// domtbl has no positive domain for.
fn classify(root: &Path, hmmer: &str, under: &[Under]) -> anyhow::Result<Vec<String>> {
    // the pairs wanted out of each unit's domtbl, which is read once through
    let mut wanted: HashMap<&str, HashMap<&str, HashMap<&str, (Vec<Domain>, u32)>>> =
        HashMap::new();
    for pair in under {
        wanted
            .entry(pair.unit.as_str())
            .or_default()
            .entry(pair.query.as_str())
            .or_default()
            .insert(pair.target.as_str(), (Vec::new(), 0));
    }

    for (unit, pairs) in &mut wanted {
        let path = manifest::dom_path(&root.join("results"), hmmer, unit);
        let file = std::fs::File::open(&path)
            .with_context(|| format!("failed to open {}", path.display()))?;
        for line in std::io::BufRead::lines(std::io::BufReader::new(file)) {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            // target, acc, tlen, query, acc, qlen, E, score, bias, #, of,
            // c-E, i-E, domain score, bias, hmm from, to, ali from, to, ...
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() < 22 {
                continue;
            }
            let Some((domains, hmm_len)) = pairs.get_mut(f[3]).and_then(|q| q.get_mut(f[0])) else {
                continue;
            };
            *hmm_len = f[5]
                .parse()
                .with_context(|| format!("bad qlen in {}", path.display()))?;
            domains.push(Domain {
                score: f[13].parse()?,
                hmm_from: f[15].parse()?,
                hmm_to: f[16].parse()?,
                ali_from: f[17].parse()?,
                ali_to: f[18].parse()?,
            });
        }
    }

    Ok(under
        .iter()
        .map(|pair| {
            let (domains, hmm_len) =
                &wanted[pair.unit.as_str()][pair.query.as_str()][pair.target.as_str()];
            match domains::measure(domains, *hmm_len) {
                Some(m) => domains::classify(&m, pair.hmmer_cutoff)
                    .0
                    .name()
                    .to_string(),
                None => "-".to_string(),
            }
        })
        .collect())
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

    /// The words of every `#= left_out` line in a written table.
    fn left_out(table: &toil::Table) -> Vec<Vec<String>> {
        table
            .meta_rows()
            .filter(|row| row.key() == "left_out")
            .map(|row| {
                (0..row.len())
                    .filter_map(|i| row.get(i))
                    .map(str::to_string)
                    .collect()
            })
            .collect()
    }

    #[test]
    fn a_family_short_of_a_cutoff_is_left_out_and_counted() {
        let dir = std::env::temp_dir().join(format!("loss-decomp-left-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        // delta has no hmmer_dom cutoff, so it is short of one
        let cutoffs = dir.join("cutoffs.tbl");
        std::fs::write(
            &cutoffs,
            "# family nail_1 mmseqs_1 hmmer_1 hmmer_dom_1\n\
             # ------ ------ -------- ------- -----------\n\
               alpha  24.5   30.0     24.5    20.0\n\
               delta  20.0   30.0     20.0    0.0\n",
        )
        .unwrap();

        // every pair seeded and none under its cutoff, so stages opens no
        // prefilter database and no domtbl. the delta rows would count as
        // truth and as kept if they were not left out
        let runs = dir.join("runs.tbl");
        std::fs::write(
            &runs,
            format!(
                "\
#= format runs 1
#= query 2 18 120
#= target 1 3 12 40
#= seed once 4.0000
#= cutoffs {} c=0
#= run hmmer hmmer 2.0000
#= run nail nail 3.0000 seeds=once
# query target           hmmer  nail   inc dom
# ----- ---------------- ------ ------ --- ---
#= shard 1
alpha MGYP000000000001 30.0   30.0   1   30.0
delta MGYP000000000002 30.0   12.0   1   30.0
delta MGYP000000000003 -      5.0    -   -
#= left_out families 1
#= left_out hmmer 1 nail 2
#= end 3
",
                cutoffs.display()
            ),
        )
        .unwrap();

        let out = dir.join("stages.tbl");
        stages(&runs, &dir, &out, 1).unwrap();

        let table = toil::Table::read(&out).unwrap();
        assert_eq!(table.rows().len(), 1);
        assert_eq!(table.get(0, "truth"), Some("1"));
        assert_eq!(table.get(0, "kept"), Some("1"));
        assert_eq!(table.get(0, "over"), Some("1"));
        assert_eq!(table.get(0, "left_out"), Some("2"));

        let want = vec![
            vec!["families".to_string(), "1".to_string()],
            ["hmmer", "1", "nail", "2"].map(str::to_string).to_vec(),
        ];
        assert_eq!(left_out(&table), want);

        let pairs = toil::Table::read(dir.join("cutoff-pairs.tbl")).unwrap();
        assert!(pairs.rows().is_empty());
        assert_eq!(left_out(&pairs), want);
    }
}
