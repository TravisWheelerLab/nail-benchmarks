//! Every score every run gave every pair, in one table.
//!
//! One row per query/target pair, one column per run, plus what hmmer scored
//! the pair and what its family's cutoffs are. This is the only thing `parse`
//! produces, and every analysis is a grouping over it: a funnel is a count
//! split by whether a pair was seeded, a sensitivity surface is a count per
//! column, a recall number is the same count per tool.
//!
//! What the columns are is read out of `manifest.tbl` rather than known here, so
//! nothing in this file can tell which pipeline it is reading.
//!
//! A run that never reported a pair holds `-` rather than a zero or a NaN.
//! Absent is not a score, and a NaN would compare false against every
//! threshold and quietly look like a real answer that missed.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{Context, bail, ensure};

use bioio::split::{self, Kind};
use bioio::tbl::{BlastTable, HitTable, HmmerTable, NailTable, hmmer::HmmerDomainTable};

use bench::manifest::{self, Manifest, Wall};
use bench::tbl;

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
            other => bail!("unknown tool {other:?} in manifest.tbl"),
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
}

/// A run and the shards it covered, which is what `collect` needs and what the
/// written table has no use for: the shards are a property of the pipeline,
/// listed once in its `#= target` lines rather than once per column.
struct Column {
    run: Run,
    /// In manifest order. Empty string for a pipeline that never named one.
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

        let manifest = Manifest::read(&dir.join("manifest.tbl"))?;

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

        let columns = runs(&manifest)?;
        ensure!(
            !columns.is_empty(),
            "no finished runs in {}/manifest.tbl",
            dir.display()
        );

        let seed_wall_s = walls(&manifest, "seed");

        let shards = shards(&columns);

        let results = dir.join("results");
        let seeded = seeds(&results, &shards)?;

        // one pass per column, keeping each column's scores while the union of
        // pairs grows
        let mut per_run: Vec<HashMap<Pair, f32>> = Vec::with_capacity(columns.len());
        let mut per_run_dom: Vec<HashMap<Pair, Vec<f32>>> = Vec::with_capacity(columns.len());
        let mut pairs: HashSet<Pair> = HashSet::new();

        for Column { run, shards } in &columns {
            let mut scores: HashMap<Pair, f32> = HashMap::new();
            let mut domains: HashMap<Pair, Vec<f32>> = HashMap::new();

            for shard in shards {
                let path = manifest::table_path(&results, &run.name, shard);
                for (pair, score) in run.tool.read(&path)? {
                    scores
                        .entry(pair)
                        .and_modify(|s| *s = s.max(score))
                        .or_insert(score);
                }

                if run.tool == Tool::Hmmer {
                    let path = manifest::dom_path(&results, &run.name, shard);
                    domains.extend(read_domains(&path)?);
                }
            }

            // which pairs are in the table is settled before any of them is
            // built: a pair earns a row by clearing a cutoff somewhere, or by
            // hmmer having reported it at all.
            //
            // hmmer's whole reported set is kept, cutoff or not, because it is
            // what the other columns are measured against -- a pair it found
            // weakly is still a pair they can be asked about. That is a rule
            // about the comparison, not about how hmmer is run: it is an
            // ordinary column everywhere else in this file.
            for (pair, score) in &scores {
                let kept = run.tool == Tool::Hmmer
                    || cutoffs.get(run.tool, &pair.0).is_some_and(|c| *score >= c);

                if kept {
                    pairs.insert(pair.clone());
                }
            }

            per_run.push(scores);
            per_run_dom.push(domains);
        }

        // sorted so the file is stable across runs and diffs mean something
        let mut pairs: Vec<Pair> = pairs.into_iter().collect();
        pairs.sort_unstable();

        let rows = pairs
            .into_iter()
            .map(|(query, target)| {
                let key = (query.clone(), target.clone());

                Row {
                    cut_nail: cutoffs.nail.get(&query).copied(),
                    cut_mmseqs: cutoffs.mmseqs.get(&query).copied(),
                    seeded: seeded.as_ref().map(|set| set.contains(&key)),
                    scores: per_run.iter().map(|r| r.get(&key).copied()).collect(),
                    domains: per_run_dom
                        .iter()
                        .map(|d| d.get(&key).cloned().unwrap_or_default())
                        .collect(),
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
            runs: columns.into_iter().map(|c| c.run).collect(),
            rows,
        })
    }

    pub fn write(&self, path: &Path) -> anyhow::Result<()> {
        // `#=` is metadata for whoever reads this back, `#` is the header a
        // person reads. stockholm draws the same line in the same place.
        let mut meta = format!(
            "#= query {} {} {}\n",
            self.query.count, self.query.residues, self.query.bytes
        );

        for (shard, size) in &self.targets {
            meta.push_str(&format!(
                "#= target {} {} {} {}\n",
                label(shard),
                size.count,
                size.residues,
                size.bytes
            ));
        }

        for (shard, wall) in &self.seed_wall_s {
            meta.push_str(&format!("#= seed {} {wall:.4}\n", label(shard)));
        }

        for run in &self.runs {
            let params: String = run
                .params
                .iter()
                .map(|(k, v)| format!(" {k}={v}"))
                .collect();
            meta.push_str(&format!(
                "#= run {} {} {:.4}{params}\n",
                run.name, run.tool, run.wall_s
            ));
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

        // a domain breakdown is as wide as the pair has domains, so those
        // columns go on the end where a ragged one costs nothing
        let dom: Vec<usize> = (0..self.runs.len())
            .filter(|&i| self.runs[i].tool == Tool::Hmmer)
            .collect();
        headers.extend(dom.iter().map(|&i| format!("{}_dom", self.runs[i].name)));

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
                row.extend(dom.iter().map(|&i| domains(&r.domains[i])));
                row
            })
            .collect();

        // the domain list is as wide as the pair has domains, so like argv in
        // manifest.tbl the last column is written as it comes and never padded
        tbl::write(
            path,
            tbl::Table {
                meta: &meta,
                headers: &headers,
                rows: &cells,
                ragged_last: true,
            },
        )
    }

    /// Reads back what [`Scores::write`] wrote.
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
        /// A tenth of the best domain, which is mgnify's threshold.
        const SIGNIFICANT: f32 = 0.1;

        let domains = &self.domains[run];
        let Some(&best) = domains.first() else {
            return 0;
        };

        match best > 0.0 {
            true => domains.iter().filter(|d| *d / best >= SIGNIFICANT).count(),
            false => domains.len(),
        }
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
        [] => dash(),
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

// ----------------------------------------------------------------- manifest.tbl

/// The columns a pipeline ran, in the order it declared them.
///
/// A run is one `name`; a name that turns up against several shards is one
/// column covering all of them, since what changed between those commands is
/// the target rather than the parameterization.
fn runs(manifest: &Manifest) -> anyhow::Result<Vec<Column>> {
    let mut out: Vec<Column> = Vec::new();
    let mut at: HashMap<String, usize> = HashMap::new();
    let mut walls: Vec<Wall> = Vec::new();

    for row in manifest.runs() {
        let name = row.get(manifest::NAME).expect("runs() filters on name");
        let shard = row.get(manifest::SHARD).unwrap_or_default().to_string();

        let i = match at.get(name) {
            Some(&i) => i,
            None => {
                let tool = row
                    .get(manifest::TOOL)
                    .with_context(|| format!("run {name:?} has no tool"))?;

                at.insert(name.to_string(), out.len());
                walls.push(Wall::default());
                out.push(Column {
                    run: Run {
                        name: name.to_string(),
                        tool: Tool::parse(tool)?,
                        wall_s: 0.0,
                        params: row.params(),
                    },
                    shards: Vec::new(),
                });
                out.len() - 1
            }
        };

        if !out[i].shards.contains(&shard) {
            out[i].shards.push(shard.clone());
        }
        walls[i].add(&shard, row);
    }

    for (column, wall) in out.iter_mut().zip(&walls) {
        column.run.wall_s = wall.total();
    }

    Ok(out)
}

/// Every shard any run covered, in the order the runs named them.
fn shards(columns: &[Column]) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for column in columns {
        for shard in &column.shards {
            if !out.contains(shard) {
                out.push(shard.clone());
            }
        }
    }
    out
}

/// What one stage of the pipeline cost, per shard.
fn walls(manifest: &Manifest, stage: &str) -> Vec<(String, f64)> {
    let mut wall = Wall::default();

    for row in manifest.stage(stage) {
        let shard = row.get(manifest::SHARD).unwrap_or_default();
        wall.add(shard, row);
    }

    wall.per_bucket()
}

fn shard_path(targets: &Path, shard: &str) -> std::path::PathBuf {
    match shard.is_empty() {
        // the one target a pipeline with no shard axis searched
        true => targets.join("target.fa"),
        false => targets.join(format!("{shard}.fa")),
    }
}

// ------------------------------------------------------------------ domains

/// The domain scores behind each of one hmmer run's hits, best first.
///
/// Only hmmer breaks a hit down this way, so this is the one thing a run of it
/// carries that a run of anything else does not. A pair in the hit table but
/// not here simply has no breakdown.
fn read_domains(path: &Path) -> anyhow::Result<HashMap<Pair, Vec<f32>>> {
    let dom = HmmerDomainTable::from_path(path, |_| true)
        .with_context(|| format!("failed to read {}", path.display()))?;

    Ok(dom
        .hits
        .into_iter()
        .map(|(pair, hit)| {
            let mut domains: Vec<f32> = hit.domains.iter().map(|d| d.score).collect();
            // best first, so the one a threshold would see is in front and the
            // rest read as the tail behind it
            domains.sort_by(|x, y| y.total_cmp(x));
            (pair, domains)
        })
        .collect())
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
