//! The first `scores.tbl` grammar: one row per pair, one column per run.
//!
//! ```text
//! # query target cut_nail cut_mmseqs seeded <run>... <run>_dom...
//! ```
//!
//! Nothing writes this any more. `parse funnel` reads it, because the funnel
//! is hit-loss's analysis and hit-loss has no table of its own yet -- it wants
//! a `seeded` column, which recall's table has no use for and does not carry.
//! Both this module and `parse funnel` go when hit-loss gets its own.

// the shape is the file's rather than one analysis's, so the parts the funnel
// does not read are still what the format has in it
#![allow(dead_code)]

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, bail, ensure};

use super::{Run, Size, Tool};

/// One query/target pair, and what everything scored it.
#[derive(Clone, Debug)]
pub struct Row {
    pub query: String,
    pub target: String,
    /// `None` for a family whose decoys gave that tool no usable cutoff.
    pub cut_nail: Option<f32>,
    pub cut_mmseqs: Option<f32>,
    /// Whether seeding found this pair. `None` when the pipeline kept no seeds.
    pub seeded: Option<bool>,
    /// One per run, in the same order.
    pub scores: Vec<Option<f32>>,
    /// One per run, in the same order: the domain scores behind that run's
    /// score, best first. Empty for every run but hmmer's, which is the only
    /// tool that breaks a hit down.
    pub domains: Vec<Vec<f32>>,
}

#[derive(Debug)]
pub struct Scores {
    pub query: Size,
    /// One per shard the runs covered, in manifest order.
    pub targets: Vec<(String, Size)>,
    /// Per shard, what seeding cost. Empty when a pipeline never seeded.
    pub seed_wall_s: Vec<(String, f64)>,
    pub runs: Vec<Run>,
    pub rows: Vec<Row>,
}

impl Scores {
    /// Reads a table written in this grammar.
    ///
    /// The `#= run` lines give the columns and the `#` header says whether a
    /// `seeded` column sits in front of them, so a row's fields are placed by
    /// what the file declares rather than by a count agreed on in advance.
    pub fn read(path: &Path) -> anyhow::Result<Scores> {
        let file =
            File::open(path).with_context(|| format!("failed to open {}", path.display()))?;

        let mut query: Option<Size> = None;
        let mut targets: Vec<(String, Size)> = Vec::new();
        let mut seed_wall_s: Vec<(String, f64)> = Vec::new();
        let mut runs: Vec<Run> = Vec::new();
        let mut seeded_col = false;
        let mut rows: Vec<Row> = Vec::new();

        for line in BufReader::new(file).lines() {
            let line = line?;

            if let Some(rest) = line.strip_prefix("#=") {
                let mut it = rest.split_whitespace();
                let Some(key) = it.next() else { continue };
                let f: Vec<&str> = it.collect();

                match key {
                    "format" => bail!(
                        "{} is a newer scores.tbl than the funnel reads",
                        path.display()
                    ),
                    "query" => query = Some(size(&f)?),
                    "target" => targets.push((shard_of(&f)?, size(&f[1..])?)),
                    "seed" => {
                        let wall = f.get(1).context("a `#= seed` line wants a wall time")?;
                        seed_wall_s.push((shard_of(&f)?, wall.parse()?));
                    }
                    "run" => runs.push(run(&f)?),
                    other => bail!("unknown `#= {other}` line in {}", path.display()),
                }
                continue;
            }

            if let Some(rest) = line.strip_prefix('#') {
                // the header names the columns; the rule under it is dashes
                let rest = rest.trim_start();
                if rest.starts_with("query ") {
                    seeded_col = rest.split_whitespace().any(|h| h == "seeded");
                }
                continue;
            }

            if line.trim().is_empty() {
                continue;
            }

            // query, target, both cutoffs, maybe seeded, one per run, then one
            // domain list per hmmer run. every cell is a single token, the
            // comma-joined domain lists included.
            let f: Vec<&str> = line.split_whitespace().collect();
            let doms = runs.iter().filter(|r| r.tool == Tool::Hmmer).count();
            let want = 4 + usize::from(seeded_col) + runs.len() + doms;
            ensure!(
                f.len() == want,
                "a row of {} has {} fields, expected {want}",
                path.display(),
                f.len()
            );

            let at = 4 + usize::from(seeded_col);
            let end = at + runs.len();

            let mut dom = f[end..].iter();
            let domains = runs
                .iter()
                .map(|r| match r.tool {
                    Tool::Hmmer => parse_domains(dom.next().expect("counted above")),
                    _ => Ok(Vec::new()),
                })
                .collect::<anyhow::Result<Vec<_>>>()?;

            rows.push(Row {
                query: f[0].to_string(),
                target: f[1].to_string(),
                cut_nail: parse_score(f[2])?,
                cut_mmseqs: parse_score(f[3])?,
                seeded: match seeded_col {
                    true => Some(f[4] == "y"),
                    false => None,
                },
                scores: f[at..end]
                    .iter()
                    .map(|s| parse_score(s))
                    .collect::<anyhow::Result<Vec<_>>>()?,
                domains,
            });
        }

        ensure!(!runs.is_empty(), "no `#= run` lines in {}", path.display());

        Ok(Scores {
            query: query.context("no `#= query` line")?,
            targets,
            seed_wall_s,
            runs,
            rows,
        })
    }

    /// Which column is hmmer's, which is what everything else is measured
    /// against.
    ///
    /// A pipeline runs one, so more than one is a table the analyses have no
    /// answer for rather than a choice to make quietly.
    pub fn hmmer(&self) -> anyhow::Result<usize> {
        let mut it = self
            .runs
            .iter()
            .enumerate()
            .filter(|(_, r)| r.tool == Tool::Hmmer);

        let (i, _) = it.next().context("no hmmer run to measure against")?;
        ensure!(it.next().is_none(), "more than one hmmer run");

        Ok(i)
    }

    /// What everything is a fraction of: the pairs hmmer found and scored over
    /// their family's cutoff.
    pub fn denominator(&self, hmmer: usize) -> usize {
        self.rows
            .iter()
            .filter(|r| r.clears(Tool::Hmmer, r.scores[hmmer]))
            .count()
    }
}

impl Row {
    /// The threshold a run of `tool` is held to on this row's family.
    pub fn cutoff(&self, tool: Tool) -> Option<f32> {
        match tool {
            // hmmer takes nail's, as it does in the calibration
            Tool::Nail | Tool::Hmmer => self.cut_nail,
            Tool::Mmseqs => self.cut_mmseqs,
        }
    }

    /// Whether a score of that tool's counts as a hit here. A pair the tool
    /// never reported does not, and neither does one whose family it has no
    /// threshold for -- an unmeasurable pair is not a found one.
    pub fn clears(&self, tool: Tool, score: Option<f32>) -> bool {
        match (self.cutoff(tool), score) {
            (Some(cutoff), Some(score)) => score >= cutoff,
            _ => false,
        }
    }

    /// How many of an hmmer column's domains carry enough of the hit to count
    /// as their own.
    ///
    /// The measure is against the best domain rather than an absolute score:
    /// what is being asked is whether the hit is one region of the sequence or
    /// several, and a weak family's several are still several.
    pub fn domain_count(&self, run: usize) -> usize {
        let domains = &self.domains[run];
        let Some(&best) = domains.first() else {
            return 0;
        };

        match best > 0.0 {
            true => domains
                .iter()
                .filter(|d| *d / best >= super::SIGNIFICANT)
                .count(),
            false => domains.len(),
        }
    }
}

fn parse_score(s: &str) -> anyhow::Result<Option<f32>> {
    match s {
        "-" => Ok(None),
        s => Ok(Some(s.parse().with_context(|| format!("bad score {s:?}"))?)),
    }
}

fn parse_domains(s: &str) -> anyhow::Result<Vec<f32>> {
    match s {
        "-" => Ok(Vec::new()),
        s => s
            .split(',')
            .map(|x| x.parse().with_context(|| format!("bad domain score {x:?}")))
            .collect(),
    }
}

/// The shard a `#= target` or `#= seed` line is about, back from its label.
fn shard_of(fields: &[&str]) -> anyhow::Result<String> {
    match *fields.first().context("a metadata line names no shard")? {
        "-" => Ok(String::new()),
        shard => Ok(shard.to_string()),
    }
}

fn size(fields: &[&str]) -> anyhow::Result<Size> {
    let [count, residues, bytes] = fields else {
        bail!("a size wants a count, a residue count and a byte count");
    };

    Ok(Size {
        count: count.parse()?,
        residues: residues.parse()?,
        bytes: bytes.parse()?,
    })
}

/// One `#= run` line: a name, a tool, a wall time, then whatever settings told
/// this run apart from the others.
fn run(fields: &[&str]) -> anyhow::Result<Run> {
    let [name, tool, wall, params @ ..] = fields else {
        bail!("a `#= run` line wants a name, a tool and a wall time");
    };

    Ok(Run {
        name: name.to_string(),
        tool: Tool::parse(tool)?,
        wall_s: wall.parse()?,
        params: params
            .iter()
            .filter_map(|p| p.split_once('='))
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect(),
    })
}
