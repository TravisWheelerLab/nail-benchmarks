//! Every score every run gave every pair, in one table.
//!
//! One row per query/target pair, one column per run, plus what hmmer scored
//! the pair and what its family's cutoffs are. This is the only thing `parse`
//! produces, and every analysis is a grouping over it: a funnel is a count
//! split by whether a pair was seeded, a sensitivity surface is a count per
//! column, a recall number is the same count per tool.
//!
//! What the columns are is read out of `runs.tbl` rather than known here, so
//! nothing in this file can tell which pipeline it is reading.
//!
//! A run that never reported a pair holds `-` rather than a zero or a NaN.
//! Absent is not a score, and a NaN would compare false against every
//! threshold and quietly look like a real answer that missed.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{Context, bail, ensure};

use bioio::split::{self, Kind};
use bioio::tbl::{BlastTable, HitTable, HmmerTable, NailTable, hmmer::HmmerDomainTable};

use crate::manifest::{self, Manifest};

type Pair = (String, String);

/// Which program produced a results table, which settles both how to read it
/// and which cutoff its scores are held against.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tool {
    Nail,
    Mmseqs,
    Hmmer,
}

impl Tool {
    fn parse(name: &str) -> anyhow::Result<Tool> {
        match name {
            "nail" => Ok(Tool::Nail),
            "mmseqs" => Ok(Tool::Mmseqs),
            "hmmer" => Ok(Tool::Hmmer),
            other => bail!("unknown tool {other:?} in runs.tbl"),
        }
    }

    /// The best score this tool gave each pair in one results table.
    ///
    /// mmseqs is read as blast: `convertalis --format-mode 0` is blast's
    /// tabular format, whatever wrote it.
    fn read(self, path: &Path) -> anyhow::Result<HashMap<Pair, f32>> {
        let tbl = match self {
            Tool::Nail => HitTable::from_path::<_, NailTable>(path),
            Tool::Mmseqs => HitTable::from_path::<_, BlastTable>(path),
            Tool::Hmmer => HitTable::from_path::<_, HmmerTable>(path),
        }
        .with_context(|| format!("failed to read {}", path.display()))?;

        let mut best: HashMap<Pair, f32> = HashMap::new();
        for hit in tbl.hits {
            // a pair can be reported more than once; the best of them is the
            // one a threshold would see
            best.entry((hit.query, hit.target))
                .and_modify(|s| *s = s.max(hit.score))
                .or_insert(hit.score);
        }

        Ok(best)
    }
}

impl fmt::Display for Tool {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Tool::Nail => "nail",
            Tool::Mmseqs => "mmseqs",
            Tool::Hmmer => "hmmer",
        };
        write!(f, "{name}")
    }
}

/// How big one side of a search is.
#[derive(Clone, Copy, Debug)]
pub struct Size {
    /// Families for the query, sequences for a target shard.
    pub count: usize,
    pub residues: u64,
    pub bytes: u64,
}

/// One column: a named run of one tool at one parameterization.
#[derive(Clone, Debug)]
pub struct Run {
    pub name: String,
    pub tool: Tool,
    /// What the column cost, summed over every shard it covered.
    pub wall_s: f64,
    /// Whatever else the commands recorded — the settings that tell this run
    /// apart from the others of the same tool.
    pub params: BTreeMap<String, String>,
    /// The shards it ran against, in manifest order. Empty string for a
    /// pipeline that never named one.
    shards: Vec<String>,
}

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
    /// One per column, in the same order.
    pub scores: Vec<Option<f32>>,
    /// hmmer's score for the whole sequence, `None` if hmmer never reported
    /// the pair. This is the one an analysis compares against; a domain score
    /// is a piece of it.
    pub hmmer_seq: Option<f32>,
    /// Every domain score hmmer gave the pair, best first.
    pub hmmer_dom: Vec<f32>,
}

#[derive(Debug)]
pub struct Scores {
    pub query: Size,
    /// One per shard the runs covered, in manifest order.
    pub targets: Vec<(String, Size)>,
    /// Per shard, what seeding and hmmer cost. Empty when a pipeline had no
    /// such stage.
    pub seed_wall_s: Vec<(String, f64)>,
    pub truth_wall_s: Vec<(String, f64)>,
    pub runs: Vec<Run>,
    pub rows: Vec<Row>,
}

/// Where a pipeline's inputs are, so `parse` can measure what was searched.
pub struct Inputs<'a> {
    pub query_hmm: &'a Path,
    pub targets: &'a Path,
}

impl Scores {
    /// Reads a finished pipeline directory into a table.
    pub fn collect(
        dir: &Path,
        inputs: Inputs<'_>,
        cutoffs_path: &Path,
        c: usize,
    ) -> anyhow::Result<Scores> {
        let cutoffs = cutoffs(cutoffs_path, c)
            .with_context(|| format!("failed to read {}", cutoffs_path.display()))?;

        let manifest = Manifest::read(&dir.join("runs.tbl"))?;

        let failed: Vec<&str> = manifest
            .failed()
            .filter_map(|row| row.get(manifest::NAME))
            .collect();
        if !failed.is_empty() {
            eprintln!(
                "warning: leaving out {} run(s) that did not finish: {}",
                failed.len(),
                failed.join(", ")
            );
        }

        let runs = runs(&manifest)?;
        ensure!(
            !runs.is_empty(),
            "no finished runs in {}/runs.tbl",
            dir.display()
        );

        // hmmer's parts run at the same time, so the longest of them is about
        // what the shard's truth cost; seeding is one command per shard
        let truth_wall_s = walls(&manifest, "truth", f64::max);
        let seed_wall_s = walls(&manifest, "seed", f64::max);

        let shards = shards(&runs);

        let results = dir.join("results");

        let truth = truth(&results, &shards)?;
        let seeded = seeds(&results, &shards)?;

        // one pass per column, keeping each column's scores while the union of
        // pairs grows
        let mut per_run: Vec<HashMap<Pair, f32>> = Vec::with_capacity(runs.len());
        let mut pairs: HashSet<Pair> = truth.keys().cloned().collect();

        for run in &runs {
            let mut scores: HashMap<Pair, f32> = HashMap::new();

            for shard in &run.shards {
                let path = manifest::table_path(&results, &run.name, shard);
                for (pair, score) in run.tool.read(&path)? {
                    scores
                        .entry(pair)
                        .and_modify(|s| *s = s.max(score))
                        .or_insert(score);
                }
            }

            // a pair earns a row by clearing the cutoff somewhere, so which
            // pairs are in the table is settled before any of them is built
            for (pair, score) in &scores {
                if cutoffs.get(run.tool, &pair.0).is_some_and(|c| *score >= c) {
                    pairs.insert(pair.clone());
                }
            }

            per_run.push(scores);
        }

        // sorted so the file is stable across runs and diffs mean something
        let mut pairs: Vec<Pair> = pairs.into_iter().collect();
        pairs.sort_unstable();

        let rows = pairs
            .into_iter()
            .map(|(query, target)| {
                let key = (query.clone(), target.clone());
                let found = truth.get(&key);

                Row {
                    cut_nail: cutoffs.nail.get(&query).copied(),
                    cut_mmseqs: cutoffs.mmseqs.get(&query).copied(),
                    seeded: seeded.as_ref().map(|set| set.contains(&key)),
                    scores: per_run.iter().map(|r| r.get(&key).copied()).collect(),
                    hmmer_seq: found.map(|(seq, _)| *seq),
                    hmmer_dom: found.map(|(_, dom)| dom.clone()).unwrap_or_default(),
                    query,
                    target,
                }
            })
            .collect();

        let targets = shards
            .iter()
            .map(|shard| {
                let size = fasta_size(&shard_path(inputs.targets, shard))?;
                Ok((shard.clone(), size))
            })
            .collect::<anyhow::Result<Vec<_>>>()?;

        Ok(Scores {
            query: hmm_size(inputs.query_hmm)?,
            targets,
            seed_wall_s,
            truth_wall_s,
            runs,
            rows,
        })
    }

    pub fn write(&self, path: &Path) -> anyhow::Result<()> {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }

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

        for (shard, size) in &self.targets {
            writeln!(
                out,
                "#= target {} {} {} {}",
                label(shard),
                size.count,
                size.residues,
                size.bytes
            )?;
        }

        for (shard, wall) in &self.seed_wall_s {
            writeln!(out, "#= seed {} {wall:.4}", label(shard))?;
        }
        for (shard, wall) in &self.truth_wall_s {
            writeln!(out, "#= truth {} {wall:.4}", label(shard))?;
        }

        for run in &self.runs {
            let params: String = run
                .params
                .iter()
                .map(|(k, v)| format!(" {k}={v}"))
                .collect();
            writeln!(
                out,
                "#= run {} {} {:.4}{params}",
                run.name, run.tool, run.wall_s
            )?;
        }

        let mut headers = vec![
            "query".to_string(),
            "target".to_string(),
            "cut_nail".to_string(),
            "cut_mmseqs".to_string(),
        ];
        let seeded = self.rows.first().is_some_and(|r| r.seeded.is_some());
        if seeded {
            headers.push("seeded".to_string());
        }
        headers.extend(self.runs.iter().map(|r| r.name.clone()));
        headers.push("hmmer_seq".to_string());
        headers.push("hmmer_dom".to_string());

        let cells: Vec<Vec<String>> = self
            .rows
            .iter()
            .map(|r| {
                let mut row = vec![
                    r.query.clone(),
                    r.target.clone(),
                    score(r.cut_nail),
                    score(r.cut_mmseqs),
                ];
                if seeded {
                    row.push(match r.seeded {
                        Some(true) => "y".to_string(),
                        Some(false) => "n".to_string(),
                        None => dash(),
                    });
                }
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
}

// ------------------------------------------------------------------- fields

fn dash() -> String {
    "-".to_string()
}

/// A shard with no name, from a pipeline that only ever searched one target.
fn label(shard: &str) -> &str {
    match shard.is_empty() {
        true => "-",
        false => shard,
    }
}

fn score(s: Option<f32>) -> String {
    match s {
        Some(s) => format!("{s:.1}"),
        None => dash(),
    }
}

/// Every domain score on one line, comma-joined so the whole list is a single
/// whitespace token however long it gets.
fn domains(scores: &[f32]) -> String {
    match scores {
        [] => dash(),
        scores => scores
            .iter()
            .map(|s| format!("{s:.1}"))
            .collect::<Vec<_>>()
            .join(","),
    }
}

// ----------------------------------------------------------------- runs.tbl

/// The columns a pipeline ran, in the order it declared them.
///
/// A run is one `name`; a name that turns up against several shards is one
/// column covering all of them, since what changed between those commands is
/// the target rather than the parameterization.
fn runs(manifest: &Manifest) -> anyhow::Result<Vec<Run>> {
    let mut out: Vec<Run> = Vec::new();
    let mut at: HashMap<String, usize> = HashMap::new();

    for row in manifest.runs() {
        let name = row.get(manifest::NAME).expect("runs() filters on name");
        let shard = row.get(manifest::SHARD).unwrap_or_default().to_string();
        let wall = row.wall_s().unwrap_or(0.0);

        match at.get(name) {
            Some(&i) => {
                let run: &mut Run = &mut out[i];
                run.wall_s += wall;
                if !run.shards.contains(&shard) {
                    run.shards.push(shard);
                }
            }
            None => {
                let tool = row
                    .get(manifest::TOOL)
                    .with_context(|| format!("run {name:?} has no tool"))?;

                at.insert(name.to_string(), out.len());
                out.push(Run {
                    name: name.to_string(),
                    tool: Tool::parse(tool)?,
                    wall_s: wall,
                    params: row.params(),
                    shards: vec![shard],
                });
            }
        }
    }

    Ok(out)
}

/// Every shard any run covered, in the order the runs named them.
fn shards(runs: &[Run]) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for run in runs {
        for shard in &run.shards {
            if !out.contains(shard) {
                out.push(shard.clone());
            }
        }
    }
    out
}

/// What one stage cost per shard, folded across the commands that shared a
/// shard — hmmer's parts run together, so `max` is what that step took.
fn walls(manifest: &Manifest, stage: &str, fold: fn(f64, f64) -> f64) -> Vec<(String, f64)> {
    let mut out: Vec<(String, f64)> = Vec::new();

    for row in manifest.stage(stage) {
        let shard = row.get(manifest::SHARD).unwrap_or_default().to_string();
        let wall = row.wall_s().unwrap_or(0.0);

        match out.iter_mut().find(|(s, _)| *s == shard) {
            Some((_, at)) => *at = fold(*at, wall),
            None => out.push((shard, wall)),
        }
    }

    out
}

fn shard_path(targets: &Path, shard: &str) -> std::path::PathBuf {
    match shard.is_empty() {
        // the one target a pipeline with no shard axis searched
        true => targets.join("target.fa"),
        false => targets.join(format!("{shard}.fa")),
    }
}

// -------------------------------------------------------------------- truth

/// What hmmer scored every pair: the whole-sequence score, and every domain
/// score behind it.
///
/// Both come off the domain table, which carries the sequence score on every
/// one of a pair's rows. The sequence table can still name a pair the domain
/// table doesn't, and such a pair belongs in the union with no score beside it.
fn truth(results: &Path, shards: &[String]) -> anyhow::Result<HashMap<Pair, (f32, Vec<f32>)>> {
    let mut out: HashMap<Pair, (f32, Vec<f32>)> = HashMap::new();

    for shard in shards {
        let tbl = manifest::truth_path(results, shard, "tbl");
        let domtbl = manifest::truth_path(results, shard, "domtbl");

        let dom = HmmerDomainTable::from_path(&domtbl, |_| true)
            .with_context(|| format!("failed to read {}", domtbl.display()))?;

        for (pair, hit) in dom.hits {
            let mut domains: Vec<f32> = hit.domains.iter().map(|d| d.score).collect();
            // best first, so the one a cutoff is held against is in front and
            // the rest read as the tail behind it
            domains.sort_by(|x, y| y.total_cmp(x));
            out.insert(pair, (hit.score, domains));
        }

        let seq = HitTable::from_path::<_, HmmerTable>(&tbl)
            .with_context(|| format!("failed to read {}", tbl.display()))?;

        for hit in seq.hits {
            out.entry((hit.query, hit.target))
                .or_insert((hit.score, Vec::new()));
        }
    }

    Ok(out)
}

/// The (query, target) pairs `--seeds-out` wrote: nail's own prf/seq column
/// order, whitespace-separated, no header. `None` when the pipeline kept none.
///
/// A sequence lives in exactly one shard, so the shards' pairs are disjoint
/// and the union of them is the whole seed set.
fn seeds(results: &Path, shards: &[String]) -> anyhow::Result<Option<HashSet<Pair>>> {
    let mut out: Option<HashSet<Pair>> = None;

    for shard in shards {
        let path = manifest::seeds_path(results, shard);
        if !path.is_file() {
            continue;
        }

        let file =
            File::open(&path).with_context(|| format!("failed to open {}", path.display()))?;
        let pairs = out.get_or_insert_with(HashSet::new);

        for line in BufReader::new(file).lines() {
            let line = line?;
            let mut fields = line.split_whitespace();
            let query = fields.next().context("a seed row has no query")?;
            let target = fields.next().context("a seed row has no target")?;
            pairs.insert((query.to_string(), target.to_string()));
        }
    }

    Ok(out)
}

// ------------------------------------------------------------------ cutoffs

pub struct Cutoffs {
    pub nail: HashMap<String, f32>,
    pub mmseqs: HashMap<String, f32>,
}

impl Cutoffs {
    /// The threshold one tool's score is held against for one family.
    ///
    /// hmmer takes nail's: nail approximates hmmer's model, and the
    /// calibration learns no threshold of hmmer's own that anything reads.
    fn get(&self, tool: Tool, family: &str) -> Option<f32> {
        match tool {
            Tool::Nail | Tool::Hmmer => self.nail.get(family).copied(),
            Tool::Mmseqs => self.mmseqs.get(family).copied(),
        }
    }
}

/// One score per family per tool, out of the decoys the calibration scored it
/// against.
///
/// A zero means the family had fewer decoys than the file has slots, so that
/// tool learned nothing about it and gets no cutoff. The two tools are kept
/// apart rather than dropped together: a family nail has a threshold for is
/// still measurable against nail, whatever mmseqs made of it.
fn cutoffs(path: &Path, c: usize) -> anyhow::Result<Cutoffs> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;

    let mut out = Cutoffs {
        nail: HashMap::new(),
        mmseqs: HashMap::new(),
    };

    for line in BufReader::new(file).lines() {
        let line = line?;

        let Some((family, rest)) = line.split_once(',') else {
            continue;
        };

        for group in rest.split("),(") {
            let group = group.trim_matches(|c| c == '(' || c == ')');
            let mut it = group.split(',');

            let Some(tool) = it.next() else { continue };
            let mut nums: Vec<f32> = it.filter_map(|x| x.parse().ok()).collect();
            // the last number is how many decoys there were, not a score
            nums.pop();

            let Some(&score) = nums.get(c) else { continue };
            if score <= 0.0 {
                continue;
            }

            match tool {
                "nail" => out.nail.insert(family.to_string(), score),
                "mmseqs" => out.mmseqs.insert(family.to_string(), score),
                // the calibration also scores hmmer, which nothing is held to
                _ => None,
            };
        }
    }

    if out.nail.is_empty() && out.mmseqs.is_empty() {
        bail!("no usable cutoffs at index {c} in {}", path.display());
    }

    Ok(out)
}

// -------------------------------------------------------------------- sizes

/// How big the query set is: models, and the positions they hold.
///
/// Counted here rather than remembered from the build, so it describes the
/// files that are present rather than the ones that were drawn.
fn hmm_size(path: &Path) -> anyhow::Result<Size> {
    let bytes = std::fs::metadata(path)
        .with_context(|| format!("failed to stat {}", path.display()))?
        .len();

    // LENG per model, which is what the index weights records by
    let models = split::index(path, Kind::Hmm)?;

    Ok(Size {
        count: models.len(),
        residues: models.iter().map(|r| r.weight).sum(),
        bytes,
    })
}

/// Records and residues in a fasta: everything that isn't a header line or
/// whitespace.
///
/// Residues is the honest unit, since neither families nor sequences are
/// uniform amounts of work.
fn fasta_size(path: &Path) -> anyhow::Result<Size> {
    let bytes = std::fs::metadata(path)
        .with_context(|| format!("failed to stat {}", path.display()))?
        .len();

    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut reader = BufReader::new(file);

    let mut buf = [0u8; 1 << 16];
    let mut count = 0usize;
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
                    count += 1;
                }
                b'\n' => in_header = false,
                _ if in_header => {}
                _ if b.is_ascii_whitespace() => {}
                _ => residues += 1,
            }
        }
    }

    Ok(Size {
        count,
        residues,
        bytes,
    })
}
