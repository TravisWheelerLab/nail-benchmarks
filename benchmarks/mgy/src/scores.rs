//! One row per query/target pair, and what the runs made of it.
//!
//! This is recall's table alone. `hit-loss` and `cloud-search` render their
//! own shapes over the same collector, so no column here has to mean
//! something for a pipeline that does not have it.
//!
//! ```text
//! #= format scores 2
//! #= query <count> <residues> <bytes>
//! #= target <shard> <count> <residues> <bytes>
//! #= cutoffs <path> c=<n>
//! #= run <name> <tool> <wall_s> [k=v ...]
//! #= pass <run name> ...
//! # query target           pass   nail   mmseqs hmmer  inc dom
//! # ----- ---------------- ------ ------ ------ ------ --- ---
//! #= shard 1
//! 2-Hacid_dh_C MGYP000522683479 NNNMMH 98.8   94.0   98.7   1   98.1
//! 2-Hacid_dh_C MGYP000715666710 nnnmmh -      -      13.7   0   9.3,2.8
//! 2-oxoacid_dh MGYP000987338150 NNNMMH 145.6  140.0  145.5  1   145.2
//! #= shard 2
//! ...
//! #= end <rows>
//! ```
//!
//! One `#= target` line per shard and one `#= run` line per run, both in
//! ledger order; `#= cutoffs` records the file and the column the pass letters
//! were judged by; `#= pass` names the run behind each character of the `pass`
//! column. A reader refuses a file that does not open `#= format scores 2`.
//!
//! Rows sit in a block per shard, sorted by (query, target) within the block.
//! A sequence lives in exactly one shard, so the blocks partition the pairs
//! and nothing has to sort the whole file at once.
//!
//! `pass` holds one character per run: `n`, `m` or `h` for the tool, uppercase
//! where that run reported the pair at or above its family's cutoff and
//! lowercase otherwise, so a pair no run passed reads `nnnmmh`. hmmer is held
//! to nail's cutoff, as it is in the calibration.
//!
//! The scores are one column per tool rather than one per run. A tool gives a
//! pair the same score wherever it reports it; what its parameterization
//! changes is which pairs it reports, and `pass` is where that is recorded.
//! `-` means no run of that tool reported the pair: absent is not a score, and
//! a zero or a NaN would compare against a threshold and look like one.
//!
//! `inc` is hmmer's tblout inclusion count and `dom` its per-domain scores in
//! domtbl order, so the k-th score is the k-th row of
//! `results/<run>.<shard>.domtbl` and the coordinates stay there.
//!
//! A pair earns a row by clearing some run's cutoff, or by hmmer having
//! reported it at all. hmmer's whole reported set is kept because it is what
//! the other tools are measured against: a pair it found weakly is still a
//! pair they can be asked about.
//!
//! Every column is padded to its width except `query`, which is unpadded
//! because a query's rows are adjacent and so line up without it, and `dom`,
//! which is as wide as the pair has domains.

pub mod read;
mod scan;
pub mod shard;
pub mod sizes;
pub mod v1;
pub mod write;

use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};

use libsail::collection::{Indexable, Iterable};
use libsail::seq::p7hmm::IndexedHmm;

use util::ledger::{self, Ledger};

pub use v1::Scores;

/// What the file opens with, and what a reader will not read past.
pub const FORMAT: &str = "#= format scores 2";

/// A tenth of the best domain, which is mgnify's threshold for a domain
/// carrying enough of a hit to count as its own.
pub const SIGNIFICANT: f32 = 0.1;

/// How many runs a pass string can hold. cloud-search's grid is the widest
/// pipeline here at 83.
pub const MAX_RUNS: usize = 128;

// ---

/// Which program produced a results table, which settles both how to read it
/// and which cutoff its scores are held against.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Tool {
    Nail,
    Mmseqs,
    Hmmer,
}

impl Tool {
    /// Every tool, in the order their score columns are written.
    pub const ALL: [Tool; 3] = [Tool::Nail, Tool::Mmseqs, Tool::Hmmer];

    pub fn parse(name: &str) -> anyhow::Result<Tool> {
        match name {
            "nail" => Ok(Tool::Nail),
            "mmseqs" => Ok(Tool::Mmseqs),
            "hmmer" => Ok(Tool::Hmmer),
            other => bail!("unknown tool {other:?} in ledger.tbl"),
        }
    }

    /// This tool's character in a `pass` string, lowercase.
    pub fn letter(self) -> u8 {
        match self {
            Tool::Nail => b'n',
            Tool::Mmseqs => b'm',
            Tool::Hmmer => b'h',
        }
    }

    /// The tool a `pass` character names, whichever case it is in.
    pub fn of_letter(letter: u8) -> Option<Tool> {
        match letter.to_ascii_lowercase() {
            b'n' => Some(Tool::Nail),
            b'm' => Some(Tool::Mmseqs),
            b'h' => Some(Tool::Hmmer),
            _ => None,
        }
    }

    /// Where this tool's score column sits among the ones a table carries.
    pub fn at(self) -> usize {
        match self {
            Tool::Nail => 0,
            Tool::Mmseqs => 1,
            Tool::Hmmer => 2,
        }
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
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
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
    /// Whatever else the commands recorded -- the settings that tell this run
    /// apart from the others of the same tool.
    pub params: BTreeMap<String, String>,
}

/// A run and the shards it covered, which is what collecting needs and what
/// the written table has no use for: the shards are a property of the
/// pipeline, listed once in its `#= target` lines rather than once per column.
pub struct Column {
    pub run: Run,
    /// In manifest order. Empty string for a pipeline that never named one.
    pub shards: Vec<String>,
}

/// The runs a pipeline declared, in the order it declared them.
///
/// The grouping and the wall clock are [`util::ledger`]'s work; what is left
/// here is which tool a column's name belongs to, since that is what says how
/// to read its table.
pub fn runs(ran: &Ledger) -> anyhow::Result<Vec<Column>> {
    let columns = ran.columns()?;

    ensure!(
        columns.len() <= MAX_RUNS,
        "{} runs, and a pass string holds {MAX_RUNS}",
        columns.len()
    );

    columns
        .into_iter()
        .map(|column| {
            // a run whose rows are all `*` says what it cost without saying
            // which targets it searched, and a table read by shard has nothing
            // to open. every row of its column would be a dash
            ensure!(
                !column.shards.is_empty(),
                "run {:?} covers no shard of its own in ledger.tbl: every row of it is `{}`",
                column.name,
                ledger::EVERY_SHARD,
            );

            Ok(Column {
                run: Run {
                    name: column.name,
                    tool: Tool::parse(&column.tool)?,
                    wall_s: column.wall_s,
                    params: column.params,
                },
                shards: column.shards,
            })
        })
        .collect()
}

/// The tools a set of runs used, in the order their columns are written.
pub fn tools(runs: &[Run]) -> Vec<Tool> {
    Tool::ALL
        .into_iter()
        .filter(|tool| runs.iter().any(|run| run.tool == *tool))
        .collect()
}

// -------------------------------------------------------------------- query

/// The query families, numbered in the order their names sort.
///
/// A family becomes a number once, here, and everything downstream carries the
/// number: it is half of a pair's sort key, and the index a cutoff is looked
/// up at. Sorting by the number is sorting by the name, which is what lets a
/// block be sorted without a string comparison in it.
pub struct Queries {
    /// Sorted, so a name's index is its rank.
    names: Vec<String>,
    at: HashMap<String, u32>,
    pub size: Size,
}

impl Queries {
    /// The families in a `query.hmm`, with what the set comes to.
    ///
    /// Counted off the file rather than remembered from the build, so it
    /// describes the models that are there.
    pub fn from_hmm(path: &Path) -> anyhow::Result<Queries> {
        let bytes = std::fs::metadata(path)
            .with_context(|| format!("failed to stat {}", path.display()))?
            .len();

        let models = IndexedHmm::open(path)
            .with_context(|| format!("failed to index {}", path.display()))?;

        let mut names: Vec<String> = Vec::with_capacity(models.len());
        let mut residues = 0u64;

        for model in models.iter() {
            names.push(model.header.name.clone());
            // LENG is the query axis of a search's matrix, and so the honest
            // measure of how much work the set is
            residues += model.header.leng as u64;
        }

        let size = Size {
            count: names.len(),
            residues,
            bytes,
        };

        names.sort_unstable();
        names.dedup();

        ensure!(
            names.len() <= Queries::MOST,
            "{} families, and a pair's key holds {}",
            names.len(),
            Queries::MOST
        );

        let at = names
            .iter()
            .enumerate()
            .map(|(i, name)| (name.clone(), i as u32))
            .collect();

        Ok(Queries { names, at, size })
    }

    /// How many families fit in the query half of a pair's key.
    pub const MOST: usize = 1 << 23;

    pub fn id(&self, name: &str) -> Option<u32> {
        self.at.get(name).copied()
    }

    pub fn name(&self, id: u32) -> &str {
        &self.names[id as usize]
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }
}

// ------------------------------------------------------------------ cutoffs

/// The score each family's hits are held to, by query id.
///
/// hmmer takes nail's: nail approximates hmmer's model, and the calibration
/// learns no threshold of hmmer's own that anything reads.
pub struct Cutoffs {
    nail: Vec<Option<f32>>,
    mmseqs: Vec<Option<f32>>,
}

impl Cutoffs {
    /// The threshold a run of `tool` is held to on one family.
    pub fn get(&self, tool: Tool, query: u32) -> Option<f32> {
        let column = match tool {
            Tool::Nail | Tool::Hmmer => &self.nail,
            Tool::Mmseqs => &self.mmseqs,
        };

        column.get(query as usize).copied().flatten()
    }

    /// One score per family per tool, out of the decoys the calibration scored
    /// it against.
    ///
    /// A zero means the family had fewer decoys than the file has slots, so
    /// that tool learned nothing about it and gets no cutoff. The two tools are
    /// kept apart rather than dropped together: a family nail has a threshold
    /// for is still measurable against nail, whatever mmseqs made of it.
    ///
    /// `cutoffs.tbl` names its columns `<tool>_1..<tool>_5` and `<tool>_n`, so
    /// the column to read is found by name rather than by counting -- which is
    /// what lets the calibration add a tool without moving anything here.
    pub fn read(path: &Path, c: usize, queries: &Queries) -> anyhow::Result<Cutoffs> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read {}", path.display()))?;

        let headers: Vec<&str> = text
            .lines()
            .find_map(|line| {
                let rest = line.strip_prefix('#')?.trim_start();
                (!rest.starts_with('-')).then(|| rest.split_whitespace().collect())
            })
            .with_context(|| format!("no header in {}", path.display()))?;

        let column = |tool: &str| -> anyhow::Result<usize> {
            let name = format!("{tool}_{}", c + 1);
            headers
                .iter()
                .position(|h| *h == name)
                .with_context(|| format!("{} has no {name} column", path.display()))
        };

        let (nail_at, mmseqs_at) = (column("nail")?, column("mmseqs")?);

        let mut out = Cutoffs {
            nail: vec![None; queries.len()],
            mmseqs: vec![None; queries.len()],
        };

        for line in text.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            let cells: Vec<&str> = line.split_whitespace().collect();
            let Some(family) = cells.first() else {
                continue;
            };

            // a family the calibration knows and this query set does not is
            // not an error: the calibration covers every family Pfam has
            let Some(id) = queries.id(family) else {
                continue;
            };

            let score = |at: usize| cells.get(at).and_then(|x| x.parse::<f32>().ok());

            // a zero is a family that tool learned nothing about
            if let Some(s) = score(nail_at).filter(|s| *s > 0.0) {
                out.nail[id as usize] = Some(s);
            }
            if let Some(s) = score(mmseqs_at).filter(|s| *s > 0.0) {
                out.mmseqs[id as usize] = Some(s);
            }
        }

        if out.nail.iter().all(Option::is_none) && out.mmseqs.iter().all(Option::is_none) {
            bail!(
                "no usable cutoffs at index {c} in {} for any of the {} families searched",
                path.display(),
                queries.len()
            );
        }

        Ok(out)
    }
}

// --------------------------------------------------------------------- meta

/// What a `scores.tbl` says about itself, above the header.
pub struct Meta {
    pub query: Size,
    /// One per shard the runs covered, in ledger order.
    pub targets: Vec<(String, Size)>,
    pub cutoffs: PathBuf,
    pub c: usize,
    pub runs: Vec<Run>,
}

impl Meta {
    pub fn write(&self, out: &mut impl std::io::Write) -> std::io::Result<()> {
        writeln!(out, "{FORMAT}")?;
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

        writeln!(
            out,
            "#= cutoffs {} c={}",
            self.cutoffs.display(),
            self.c
        )?;

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

        let names: Vec<&str> = self.runs.iter().map(|run| run.name.as_str()).collect();
        writeln!(out, "#= pass {}", names.join(" "))
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
            .filter(|(_, run)| run.tool == Tool::Hmmer);

        let (i, _) = it.next().context("no hmmer run to measure against")?;
        ensure!(it.next().is_none(), "more than one hmmer run");

        Ok(i)
    }

    /// The tools whose score columns the table carries.
    pub fn tools(&self) -> Vec<Tool> {
        tools(&self.runs)
    }
}

/// A [`Meta`] as its lines arrive, since a reader meets them one at a time.
#[derive(Default)]
pub struct Preamble {
    query: Option<Size>,
    targets: Vec<(String, Size)>,
    cutoffs: Option<PathBuf>,
    c: Option<usize>,
    runs: Vec<Run>,
    pass: Vec<String>,
}

impl Preamble {
    /// Take in one `#=` line's key and fields.
    pub fn absorb(&mut self, key: &str, fields: &[&str]) -> anyhow::Result<()> {
        match key {
            "query" => self.query = Some(size(fields)?),
            "target" => self.targets.push((shard_of(fields)?, size(&fields[1..])?)),
            "cutoffs" => {
                let [path, rest @ ..] = fields else {
                    bail!("a `#= cutoffs` line wants a path");
                };

                self.cutoffs = Some(PathBuf::from(path));
                self.c = rest
                    .iter()
                    .find_map(|field| field.strip_prefix("c="))
                    .map(str::parse)
                    .transpose()?;
            }
            "run" => self.runs.push(run(fields)?),
            "pass" => self.pass = fields.iter().map(|name| name.to_string()).collect(),
            other => bail!("unknown `#= {other}` line"),
        }

        Ok(())
    }

    /// The preamble, once the header has been reached.
    pub fn finish(self) -> anyhow::Result<Meta> {
        let meta = Meta {
            query: self.query.context("no `#= query` line")?,
            targets: self.targets,
            cutoffs: self.cutoffs.context("no `#= cutoffs` line")?,
            c: self.c.context("no `c=` on the `#= cutoffs` line")?,
            runs: self.runs,
        };

        ensure!(!meta.runs.is_empty(), "no `#= run` lines");
        ensure!(!meta.targets.is_empty(), "no `#= target` lines");

        // the pass string is read by position, so a legend that disagrees with
        // the runs is a file that cannot be read rather than one to guess at
        let names: Vec<&str> = meta.runs.iter().map(|run| run.name.as_str()).collect();
        ensure!(
            self.pass == names,
            "`#= pass` names {:?}, the runs are {:?}",
            self.pass,
            names
        );

        Ok(meta)
    }
}

/// A shard with no name, from a pipeline that only ever searched one target.
pub fn label(shard: &str) -> &str {
    match shard.is_empty() {
        true => "-",
        false => shard,
    }
}

/// The shard a `#= target` line is about, back from its label.
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
