//! Every score every parameterization gave every pair.
//!
//! One row per query/target pair, one column per cell on the grid, plus what
//! hmmer scored the pair and what its family's cutoff is. The rows are the
//! union of everything anything reported, so a pair is here if hmmer found it
//! or if any single cell did.
//!
//! This is the expensive half of the analysis and the half unlikely to change
//! shape: reading it means reading every results table once. grid.tbl is a
//! function of it, so a different statistic or a different table shape is a
//! re-parse rather than a re-run.
//!
//! A cell that never reported a pair holds `None` rather than a zero or a NaN.
//! Absent is not a score, and a NaN would compare false against every threshold
//! and quietly look like a real answer that missed.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{Context, bail, ensure};

use bioio::split::{self, Kind};
use bioio::tbl::{HitTable, HmmerTable, NailTable, hmmer::HmmerDomainTable};

use crate::cell::Cell;

/// How big one side of the search is.
#[derive(Clone, Copy, Debug)]
pub struct Size {
    /// Families for the query, sequences for the target.
    pub count: usize,
    pub residues: u64,
    pub bytes: u64,
}

/// One cell on the grid, and what it cost.
#[derive(Clone, Copy, Debug)]
pub struct Column {
    pub cell: Cell,
    pub wall_s: f64,
}

/// One query/target pair, and what everything scored it.
///
/// hmmer's two scores sit at the end because the domain list is as long as the
/// pair has domains, and a column that changes width row to row can only be the
/// last one.
#[derive(Clone, Debug)]
pub struct Row {
    pub query: String,
    pub target: String,
    /// `None` for a family with no usable cutoff, which can't count either way.
    pub cutoff: Option<f32>,
    /// One per column, in the same order.
    pub scores: Vec<Option<f32>>,
    /// hmmer's score for the whole sequence, `None` if hmmer never reported the
    /// pair. This is the one the analysis compares against; see [`Row::hmmer`].
    pub hmmer_seq: Option<f32>,
    /// Every domain score hmmer gave the pair, best first. Empty if hmmer never
    /// reported it. Carried so the comparison can be tried the other way
    /// without re-reading every results table.
    pub hmmer_dom: Vec<f32>,
}

impl Row {
    /// The hmmer score a cutoff is held against: the one for the whole
    /// sequence.
    ///
    /// This is hmmer's actual answer about the pair -- what it reports, ranks
    /// by, and computes an E-value from. A domain score is a piece of it.
    ///
    /// The two don't reconcile and no setting makes them. Domain scores sum to
    /// less than the sequence score by an amount that grows with sequence
    /// length -- around 0.25 bits under 100 residues and 1.75 over 800 -- which
    /// is the flanking cost the whole-sequence score carries and an
    /// envelope-restricted domain score doesn't. Opening --domT only makes it
    /// worse, by adding weak domains that overshoot.
    pub fn hmmer(&self) -> Option<f32> {
        self.hmmer_seq
    }

    /// Whether a score clears this pair's cutoff. A pair whose family has no
    /// cutoff clears nothing, since there is no threshold to hold it to.
    pub fn clears(&self, score: Option<f32>) -> bool {
        match (self.cutoff, score) {
            (Some(cutoff), Some(score)) => score >= cutoff,
            _ => false,
        }
    }
}

#[derive(Debug)]
pub struct Scores {
    pub query: Size,
    pub target: Size,
    pub seed_wall_s: f64,
    pub hmmer_wall_s: f64,
    pub columns: Vec<Column>,
    pub rows: Vec<Row>,
}

impl Scores {
    /// Reads a finished benchmark directory into a table.
    pub fn collect(bench: &Path, cutoffs_path: &Path, c: usize) -> anyhow::Result<Scores> {
        let cutoffs = cutoffs(cutoffs_path, c)
            .with_context(|| format!("failed to read {}", cutoffs_path.display()))?;

        let runs = rows(&bench.join("runs.tbl"))?;
        let columns = columns(&runs)?;
        ensure!(
            !columns.is_empty(),
            "no cells in runs.tbl; has `cloud-search run` finished?"
        );

        // hmmer's parts run at the same time, so the longest of them is about
        // what the step took. its own summary row can't be used: a step row
        // carries no fields, so there is no stage on it to find it by.
        let hmmer_wall_s = walls(&runs, "truth").into_iter().fold(0.0, f64::max);
        let seed_wall_s = walls(&runs, "seed").into_iter().fold(0.0, f64::max);

        // the sequence score and every domain score, both off the domain table
        // -- it carries the sequence score on every one of a pair's rows
        let truth = bench.join("truth");
        let mut hmmer: HashMap<(String, String), (f32, Vec<f32>)> = HashMap::new();
        let dom = HmmerDomainTable::from_path(truth.join("hmmer.domtbl"), |_| true)?;
        for (pair, hit) in dom.hits {
            let mut domains: Vec<f32> = hit.domains.iter().map(|d| d.score).collect();
            // best first, so the one a cutoff is held against is the one in
            // front and the rest read as the tail behind it
            domains.sort_by(|x, y| y.total_cmp(x));
            hmmer.insert(pair, (hit.score, domains));
        }

        // the sequence table can name a pair the domain table doesn't, and it
        // still belongs in the union even with no score to put beside it
        let seq = HitTable::from_path::<_, HmmerTable>(truth.join("hmmer.tbl"))?;
        let mut pairs: HashSet<(String, String)> =
            seq.hits.into_iter().map(|h| (h.query, h.target)).collect();
        pairs.extend(hmmer.keys().cloned());

        // one pass per cell, keeping each cell's scores while the union grows
        let results = bench.join("results");
        let mut per_cell: Vec<HashMap<(String, String), f32>> = Vec::with_capacity(columns.len());

        for column in &columns {
            let path = results.join(format!("{}.tbl", column.cell.label()));
            let tbl = HitTable::from_path::<_, NailTable>(&path)
                .with_context(|| format!("failed to read {}", path.display()))?;

            let mut scores: HashMap<(String, String), f32> = HashMap::new();
            for hit in tbl.hits {
                let pair = (hit.query, hit.target);
                // a pair can be reported more than once; the best of them is
                // the one a threshold would see
                scores
                    .entry(pair.clone())
                    .and_modify(|s| *s = s.max(hit.score))
                    .or_insert(hit.score);
                pairs.insert(pair);
            }

            per_cell.push(scores);
        }

        // sorted so the file is stable across runs and diffs mean something
        let mut pairs: Vec<(String, String)> = pairs.into_iter().collect();
        pairs.sort_unstable();

        let rows = pairs
            .into_iter()
            .map(|(query, target)| {
                let key = (query.clone(), target.clone());
                let found = hmmer.get(&key);

                Row {
                    cutoff: cutoffs.get(&query).copied(),
                    scores: per_cell.iter().map(|c| c.get(&key).copied()).collect(),
                    hmmer_seq: found.map(|(seq, _)| *seq),
                    hmmer_dom: found.map(|(_, dom)| dom.clone()).unwrap_or_default(),
                    query,
                    target,
                }
            })
            .collect();

        Ok(Scores {
            query: sizes(&bench.join("queries/query.hmm"), Axis::Query)?,
            target: sizes(&bench.join("targets/target.fa"), Axis::Target)?,
            seed_wall_s,
            hmmer_wall_s,
            columns,
            rows,
        })
    }

    pub fn write(&self, path: &Path) -> anyhow::Result<()> {
        let file =
            File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
        let mut out = BufWriter::new(file);

        // `#=` is metadata for whoever reads this back, `#` is the header a
        // person reads. stockholm draws the same line in the same place.
        writeln!(
            out,
            "#= query {} {} {}",
            self.query.count, self.query.residues, self.query.bytes
        )?;
        writeln!(
            out,
            "#= target {} {} {}",
            self.target.count, self.target.residues, self.target.bytes
        )?;
        writeln!(out, "#= seed_wall_s {:.4}", self.seed_wall_s)?;
        writeln!(out, "#= hmmer_wall_s {:.4}", self.hmmer_wall_s)?;

        for column in &self.columns {
            let (a, b) = match column.cell {
                Cell::Pruned { a, b } => (format!("{a:.1}"), format!("{b:.1}")),
                Cell::Full => ("-".to_string(), "-".to_string()),
            };
            writeln!(
                out,
                "#= cell {} {a} {b} {:.4}",
                column.cell.label(),
                column.wall_s
            )?;
        }

        let mut headers = vec![
            "query".to_string(),
            "target".to_string(),
            "cutoff".to_string(),
        ];
        headers.extend(self.columns.iter().map(|c| c.cell.label()));
        headers.push("hmmer_seq".to_string());
        headers.push("hmmer_dom".to_string());

        let cells: Vec<Vec<String>> = self
            .rows
            .iter()
            .map(|r| {
                let mut row = vec![r.query.clone(), r.target.clone(), score(r.cutoff)];
                row.extend(r.scores.iter().map(|s| score(*s)));
                row.push(score(r.hmmer_seq));
                row.push(domains(&r.hmmer_dom));
                row
            })
            .collect();

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

        // the domain list is as wide as the pair has domains, so like argv in
        // runs.tbl the last column is written as it comes and never padded
        let last = headers.len() - 1;

        write!(out, "#")?;
        for (i, (h, &w)) in headers.iter().zip(&widths).enumerate() {
            match i == last {
                true => write!(out, " {h}")?,
                false => write!(out, " {h:<w$}")?,
            }
        }
        write!(out, "\n#")?;
        for (i, &w) in widths.iter().enumerate() {
            let w = match i == last {
                true => headers[last].len(),
                false => w,
            };
            write!(out, " {}", "-".repeat(w))?;
        }
        writeln!(out)?;

        for row in &cells {
            // the two leading spaces the `# ` takes on a header line, so the
            // columns sit under their names rather than beside them
            write!(out, "  ")?;
            for (i, (c, &w)) in row.iter().zip(&widths).enumerate() {
                if i > 0 {
                    write!(out, " ")?;
                }
                match i == last {
                    true => write!(out, "{c}")?,
                    false => write!(out, "{c:<w$}")?,
                }
            }
            writeln!(out)?;
        }

        out.flush()?;
        Ok(())
    }

    pub fn read(path: &Path) -> anyhow::Result<Scores> {
        let file =
            File::open(path).with_context(|| format!("failed to open {}", path.display()))?;

        let mut meta: BTreeMap<String, Vec<String>> = BTreeMap::new();
        let mut columns: Vec<Column> = Vec::new();
        let mut rows: Vec<Row> = Vec::new();

        for line in BufReader::new(file).lines() {
            let line = line?;

            if let Some(rest) = line.strip_prefix("#=") {
                let mut it = rest.split_whitespace();
                let Some(key) = it.next() else { continue };
                let rest: Vec<String> = it.map(str::to_string).collect();

                match key {
                    "cell" => columns.push(column(&rest)?),
                    _ => {
                        meta.insert(key.to_string(), rest);
                    }
                }
                continue;
            }

            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            // query, target, cutoff, one per cell, then hmmer's two. the domain
            // list is comma-joined, so it is one token however many it holds.
            let f: Vec<&str> = line.split_whitespace().collect();
            let want = columns.len() + 5;
            ensure!(
                f.len() == want,
                "a row has {} fields, expected {want}",
                f.len()
            );

            let end = columns.len() + 3;

            rows.push(Row {
                query: f[0].to_string(),
                target: f[1].to_string(),
                cutoff: parse_score(f[2])?,
                scores: f[3..end]
                    .iter()
                    .map(|s| parse_score(s))
                    .collect::<anyhow::Result<Vec<_>>>()?,
                hmmer_seq: parse_score(f[end])?,
                hmmer_dom: parse_domains(f[end + 1])?,
            });
        }

        ensure!(!columns.is_empty(), "no cells in {}", path.display());

        Ok(Scores {
            query: size(&meta, "query")?,
            target: size(&meta, "target")?,
            seed_wall_s: number(&meta, "seed_wall_s")?,
            hmmer_wall_s: number(&meta, "hmmer_wall_s")?,
            columns,
            rows,
        })
    }
}

// ------------------------------------------------------------------- fields

fn score(s: Option<f32>) -> String {
    match s {
        Some(s) => format!("{s:.1}"),
        None => "-".to_string(),
    }
}

fn parse_score(s: &str) -> anyhow::Result<Option<f32>> {
    match s {
        "-" => Ok(None),
        s => Ok(Some(s.parse().with_context(|| format!("bad score {s:?}"))?)),
    }
}

/// Every domain score on one line, comma-joined so the whole list is a single
/// whitespace token however long it gets.
fn domains(scores: &[f32]) -> String {
    match scores {
        [] => "-".to_string(),
        scores => scores
            .iter()
            .map(|s| format!("{s:.1}"))
            .collect::<Vec<_>>()
            .join(","),
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

fn column(fields: &[String]) -> anyhow::Result<Column> {
    let [_label, a, b, wall] = fields else {
        bail!("a `#= cell` line wants a label, an A, a B and a wall time");
    };

    let cell = match (a.as_str(), b.as_str()) {
        ("-", "-") => Cell::Full,
        (a, b) => Cell::Pruned {
            a: a.parse()?,
            b: b.parse()?,
        },
    };

    Ok(Column {
        cell,
        wall_s: wall.parse()?,
    })
}

fn size(meta: &BTreeMap<String, Vec<String>>, key: &str) -> anyhow::Result<Size> {
    let v = meta
        .get(key)
        .with_context(|| format!("no `#= {key}` line"))?;

    let [count, residues, bytes] = v.as_slice() else {
        bail!("`#= {key}` wants a count, a residue count and a byte count");
    };

    Ok(Size {
        count: count.parse()?,
        residues: residues.parse()?,
        bytes: bytes.parse()?,
    })
}

fn number(meta: &BTreeMap<String, Vec<String>>, key: &str) -> anyhow::Result<f64> {
    meta.get(key)
        .and_then(|v| v.first())
        .with_context(|| format!("no `#= {key}` line"))?
        .parse()
        .with_context(|| format!("`#= {key}` isn't a number"))
}

// ------------------------------------------------------------------ cutoffs

type Cutoffs = HashMap<String, f32>;

/// One score per family, out of the decoys mgnify scored it against.
///
/// A family is kept only when both nail and mmseqs got a nonzero cutoff, which
/// is mgnify's rule and is copied rather than improved so the numbers stay
/// comparable to its own. Nothing here runs mmseqs, and a zero means the family
/// had fewer decoys than the file has slots -- so this quietly drops the
/// cleanest families. It costs 36 of Pfam's 20795.
fn cutoffs(path: &Path, c: usize) -> anyhow::Result<Cutoffs> {
    let reader = BufReader::new(File::open(path)?);
    let mut out = Cutoffs::new();

    for line in reader.lines() {
        let line = line?;

        let Some((family, rest)) = line.split_once(',') else {
            continue;
        };

        let groups: Vec<(&str, Vec<f32>)> = rest
            .split("),(")
            .map(|g| g.trim_matches(|c| c == '(' || c == ')'))
            .map(|g| {
                let mut it = g.split(',');
                let tool = it.next().unwrap_or_default();
                let mut nums: Vec<f32> = it.filter_map(|x| x.parse().ok()).collect();
                // the last number is how many decoys there were, not a score
                nums.pop();
                (tool, nums)
            })
            .collect();

        let nail = groups.iter().find(|(t, _)| *t == "nail");
        let mmseqs = groups.iter().find(|(t, _)| *t == "mmseqs");

        let (Some((_, nail)), Some((_, mmseqs))) = (nail, mmseqs) else {
            continue;
        };

        if let (Some(&n), Some(&m)) = (nail.get(c), mmseqs.get(c))
            && n > 0.0
            && m > 0.0
        {
            out.insert(family.to_string(), n);
        }
    }

    if out.is_empty() {
        bail!("no usable cutoffs at index {c}");
    }

    Ok(out)
}

// ----------------------------------------------------------------- runs.tbl

/// The rows a sweep wrote, keyed by the header's column names.
///
/// Splitting on whitespace works for every kind of row, including the
/// continuation rows a batched step lists its commands on: `||` sits in the
/// step cell and takes one token, so the fields behind it still line up. Rows
/// with no `stage` are step summaries, and the callers filter them out by
/// asking for a stage.
///
/// `argv` is last and full of spaces, so it comes back as its first word. It
/// is not read.
fn rows(path: &Path) -> anyhow::Result<Vec<HashMap<String, String>>> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;

    let mut names: Vec<String> = Vec::new();
    let mut out = Vec::new();

    for line in BufReader::new(file).lines() {
        let line = line?;

        if let Some(rest) = line.strip_prefix('#') {
            // the first comment is the header, the second is the rule under it
            if names.is_empty() && !rest.trim_start().starts_with('-') {
                names = rest.split_whitespace().map(str::to_string).collect();
            }
            continue;
        }

        if line.trim().is_empty() {
            continue;
        }

        out.push(
            names
                .iter()
                .cloned()
                .zip(line.split_whitespace().map(str::to_string))
                .collect(),
        );
    }

    ensure!(!names.is_empty(), "no header in {}", path.display());
    Ok(out)
}

fn field<'a>(row: &'a HashMap<String, String>, key: &str) -> Option<&'a str> {
    row.get(key).map(String::as_str)
}

fn walls(rows: &[HashMap<String, String>], stage: &str) -> Vec<f64> {
    rows.iter()
        .filter(|r| field(r, "stage") == Some(stage))
        .filter_map(|r| field(r, "wall(s)")?.parse().ok())
        .collect()
}

/// The grid as the sweep actually ran it, in the order it ran.
fn columns(rows: &[HashMap<String, String>]) -> anyhow::Result<Vec<Column>> {
    rows.iter()
        .filter(|r| field(r, "stage") == Some("cell"))
        .map(|r| {
            let cell = match (field(r, "A"), field(r, "B")) {
                (Some("-"), Some("-")) => Cell::Full,
                (Some(a), Some(b)) => Cell::Pruned {
                    a: a.parse()?,
                    b: b.parse()?,
                },
                _ => bail!("a cell row has no A and B"),
            };

            Ok(Column {
                cell,
                wall_s: field(r, "wall(s)")
                    .context("a cell row has no wall time")?
                    .parse()?,
            })
        })
        .collect()
}

// -------------------------------------------------------------------- sizes

enum Axis {
    Query,
    Target,
}

/// How big one side of the search actually is.
///
/// Counted here rather than remembered from the build, so it describes the
/// files that are present rather than the ones that were drawn. Reading all of
/// Pfam costs well under a second, which is nothing next to the sweep it is
/// describing.
///
/// Residues is the honest unit: neither families nor sequences are uniform
/// amounts of work, so it is what makes a runtime here comparable to a runtime
/// from a differently sized benchmark.
fn sizes(path: &Path, axis: Axis) -> anyhow::Result<Size> {
    let bytes = std::fs::metadata(path)
        .with_context(|| format!("failed to stat {}", path.display()))?
        .len();

    let (count, residues) = match axis {
        // LENG per model, which is what the index weights records by
        Axis::Query => {
            let models = split::index(path, Kind::Hmm)?;
            (models.len(), models.iter().map(|r| r.weight).sum())
        }
        Axis::Target => fasta(path)?,
    };

    Ok(Size {
        count,
        residues,
        bytes,
    })
}

/// Records and residues in a fasta: everything that isn't a header line or
/// whitespace.
fn fasta(path: &Path) -> anyhow::Result<(usize, u64)> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut reader = BufReader::new(file);

    let mut buf = [0u8; 1 << 16];
    let mut records = 0usize;
    let mut residues = 0u64;
    let mut in_header = false;

    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }

        for &b in &buf[..n] {
            match b {
                b'>' => {
                    in_header = true;
                    records += 1;
                }
                b'\n' => in_header = false,
                _ if in_header => {}
                _ if b.is_ascii_whitespace() => {}
                _ => residues += 1,
            }
        }
    }

    Ok((records, residues))
}
