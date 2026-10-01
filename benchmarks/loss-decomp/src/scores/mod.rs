//! One row per query/target pair, and what the runs made of it.
//!
//! [`runs`] writes and reads the table, a score per run. What is in this
//! module is what it does not own: the tools, the runs a ledger declares, the
//! query numbering, the cutoffs, and the preamble. [`shard`] reads a shard's
//! results into pairs and [`collect`] turns shards into a file; [`frame`]
//! reads one back.
//!
//! A pair earns a row by clearing some run's cutoff. The most sensitive run
//! is the ceiling, and every other run is measured against what it kept.

pub mod analyze;
pub mod collect;
pub mod depth;
pub mod frame;
pub mod parse;
pub mod runs;
mod scan;
pub mod shard;

use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};

use libsail::format::Format;
use libsail::index::Reader;
use libsail::seq::p7hmm::leng_of;

use util::ledger::{self, Ledger};
use util::set::Set;

// ---

/// Which program produced a results table, which settles both how to read it
/// and which cutoff its scores are held against.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Tool {
    Nail,
    Mmseqs,
}

impl Tool {
    pub fn parse(name: &str) -> anyhow::Result<Tool> {
        match name {
            "nail" => Ok(Tool::Nail),
            "mmseqs" => Ok(Tool::Mmseqs),
            other => bail!("unknown tool {other:?} in ledger.tbl"),
        }
    }

    /// This tool's character in a `pass` string, lowercase.
    pub fn letter(self) -> u8 {
        match self {
            Tool::Nail => b'n',
            Tool::Mmseqs => b'm',
        }
    }
}

impl fmt::Display for Tool {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Tool::Nail => "nail",
            Tool::Mmseqs => "mmseqs",
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

/// What the set says each of these units came to.
///
/// Counting a shard means reading it, and at a thousand shards of three
/// gigabytes that is the whole benchmark read a second time to fill in a
/// metadata line. The build counted as it dealt and wrote the counts into the
/// manifest, so this is a lookup.
pub fn target_sizes(set: &Set, shards: &[String]) -> anyhow::Result<Vec<(String, Size)>> {
    shards
        .iter()
        .map(|shard| {
            let unit = set
                .units()
                .find(|u| u.name() == shard)
                .with_context(|| format!("the set has no unit {shard:?}"))?;

            Ok((
                shard.clone(),
                Size {
                    count: unit.number("seqs")? as usize,
                    residues: unit.number("residues")?,
                    bytes: unit.number("bytes")?,
                },
            ))
        })
        .collect()
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
    /// Which seed list this run replayed, `None` where it did not replay one.
    ///
    /// Held apart from `params` because it says where to read something rather
    /// than what was swept: as a setting it would become a column of every
    /// summary and a term in every figure's parameter space, and it is neither.
    pub seeds: Option<String>,
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

            let mut params = column.params;
            let seeds = params.remove(util::manifest::SEEDS);

            Ok(Column {
                run: Run {
                    name: column.name,
                    tool: Tool::parse(&column.tool)?,
                    wall_s: column.wall_s,
                    params,
                    seeds,
                },
                shards: column.shards,
            })
        })
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

        // one streaming pass over the framed bytes: a name and a LENG are
        // both read off a model's text, and parsing 20,795 of them to reach
        // two fields allocates every match, insert and transition row on the
        // way past
        let file = std::fs::File::open(path)
            .with_context(|| format!("failed to open {}", path.display()))?;
        let mut models = Reader::new(std::io::BufReader::new(file), Format::Hmm);

        let mut names: Vec<String> = Vec::new();
        let mut residues = 0u64;

        while models
            .advance()
            .with_context(|| format!("failed to read {}", path.display()))?
        {
            let model = models.record();

            let name = libsail::seq::name_of(Format::Hmm, model)
                .with_context(|| format!("a model in {} has no NAME", path.display()))?;
            names.push(String::from_utf8_lossy(name).into_owned());

            // LENG is the query axis of a search's matrix, and so the honest
            // measure of how much work the set is
            residues += leng_of(model)
                .with_context(|| format!("a model in {} has no LENG", path.display()))?
                as u64;
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
pub struct Cutoffs {
    nail: Vec<Option<f32>>,
    mmseqs: Vec<Option<f32>>,
}

impl Cutoffs {
    /// The threshold a run of `tool` is held to on one family.
    pub fn get(&self, tool: Tool, query: u32) -> Option<f32> {
        let column = match tool {
            Tool::Nail => &self.nail,
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
        let mut out = Cutoffs {
            nail: vec![None; queries.len()],
            mmseqs: vec![None; queries.len()],
        };

        parse(path, c, |family, nail, mmseqs| {
            // a family the calibration knows and this query set does not is
            // not an error: the calibration covers every family Pfam has
            if let Some(id) = queries.id(family) {
                out.nail[id as usize] = nail;
                out.mmseqs[id as usize] = mmseqs;
            }
        })?;

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

/// The same table keyed by the family name a row carries.
///
/// A reader has the names the table was written with and not the query set
/// the search numbered them by, so it cannot index [`Cutoffs`].
pub struct Named(HashMap<String, (Option<f32>, Option<f32>)>);

impl Named {
    /// Both tools' thresholds on one family, for a reader that asks once per
    /// family and then answers per run.
    pub fn pair(&self, family: &str) -> (Option<f32>, Option<f32>) {
        self.0.get(family).copied().unwrap_or((None, None))
    }

    pub fn read(path: &Path, c: usize) -> anyhow::Result<Named> {
        let mut out = HashMap::new();

        parse(path, c, |family, nail, mmseqs| {
            out.insert(family.to_string(), (nail, mmseqs));
        })?;

        ensure!(
            out.values()
                .any(|(nail, mmseqs)| nail.is_some() || mmseqs.is_some()),
            "no usable cutoffs at index {c} in {}",
            path.display()
        );

        Ok(Named(out))
    }
}

/// Read `cutoffs.tbl`, handing each family its two thresholds at index `c`.
///
/// `cutoffs.tbl` names its columns `<tool>_1..<tool>_5` and `<tool>_n`, so the
/// column to read is found by name rather than by counting -- which is what
/// lets the calibration add a tool without moving anything here.
fn parse(
    path: &Path,
    c: usize,
    mut row: impl FnMut(&str, Option<f32>, Option<f32>),
) -> anyhow::Result<()> {
    let table =
        toil::Table::read(path).with_context(|| format!("failed to read {}", path.display()))?;

    let column = |tool: &str| -> anyhow::Result<usize> {
        let name = format!("{tool}_{}", c + 1);
        table
            .index(&name)
            .with_context(|| format!("{} has no {name} column", path.display()))
    };

    let (nail_at, mmseqs_at) = (column("nail")?, column("mmseqs")?);

    for cells in table.rows() {
        let Some(family) = cells.get(0) else {
            continue;
        };

        // a zero is a family that tool learned nothing about
        let score = |at: usize| {
            cells
                .get(at)
                .and_then(|x| x.parse::<f32>().ok())
                .filter(|s| *s > 0.0)
        };

        row(family, score(nail_at), score(mmseqs_at));
    }

    Ok(())
}

// --------------------------------------------------------------------- meta

/// What a `scores.tbl` says about itself, above the header.
pub struct Meta {
    pub query: Size,
    /// One per shard the runs covered, in ledger order.
    pub targets: Vec<(String, Size)>,
    /// Per shard, what seeding cost. Empty where a pipeline never seeded, so
    /// that seeding taking no time and there being no seeding stay different
    /// answers.
    pub seeds: Vec<(String, f64)>,
    pub cutoffs: PathBuf,
    pub c: usize,
    pub runs: Vec<Run>,
    /// The binaries that produced the results this was read from, carried
    /// through from the run's ledger. Empty for a table written before they
    /// were recorded.
    pub tools: Vec<util::tools::Identity>,
}

impl Meta {
    /// `format` is the table's own `#= format` line. Everything under it is
    /// the same whatever the columns turn out to be, which is why the two
    /// tables share a preamble and not a schema.
    pub fn write<W: std::io::Write>(
        &self,
        format: &str,
        out: &mut toil::Stream<W>,
    ) -> std::io::Result<()> {
        out.meta("format", [format])?;

        let size = |size: &Size| {
            [
                toil::Cell::from(size.count),
                size.residues.into(),
                size.bytes.into(),
            ]
        };
        out.meta("query", size(&self.query))?;

        for (shard, size_of) in &self.targets {
            let words = [toil::Cell::from(shard)].into_iter().chain(size(size_of));
            out.meta("target", words)?;
        }

        for id in &self.tools {
            out.meta("tool", [&id.name, &id.version, &id.hash])?;
        }

        for (name, wall) in &self.seeds {
            out.meta("seed", [name.clone(), format!("{wall:.4}")])?;
        }

        out.meta(
            "cutoffs",
            [self.cutoffs.display().to_string(), format!("c={}", self.c)],
        )?;

        for run in &self.runs {
            // written with the settings and read back out of them, so the
            // line stays one shape and a reader needs no new field
            let words = [
                run.name.clone(),
                run.tool.to_string(),
                format!("{:.4}", run.wall_s),
            ]
            .into_iter()
            .chain(
                run.seeds
                    .iter()
                    .map(|name| format!("{}={name}", util::manifest::SEEDS)),
            )
            .chain(run.params.iter().map(|(k, v)| format!("{k}={v}")));

            out.meta("run", words)?;
        }

        Ok(())
    }

    /// Which column is the ceiling, which every other run is measured
    /// against: the run seeded at the highest sensitivity.
    pub fn ceiling(&self) -> anyhow::Result<usize> {
        let mut at: Vec<(usize, f64)> = Vec::with_capacity(self.runs.len());
        for (i, run) in self.runs.iter().enumerate() {
            let s = run
                .params
                .get("s")
                .with_context(|| format!("run {:?} records no `s` setting", run.name))?;
            let s: f64 = s
                .parse()
                .with_context(|| format!("run {:?} has s={s:?}", run.name))?;
            at.push((i, s));
        }

        let top = at.iter().map(|&(_, s)| s).fold(f64::NEG_INFINITY, f64::max);
        let mut it = at.iter().filter(|&&(_, s)| s == top);

        let &(i, _) = it.next().context("no runs to take a ceiling from")?;
        ensure!(
            it.next().is_none(),
            "two runs seeded at s={top}; no single run is the ceiling"
        );

        Ok(i)
    }
}

/// A [`Meta`] as its lines arrive, since a reader meets them one at a time.
#[derive(Default)]
pub struct Preamble {
    query: Option<Size>,
    targets: Vec<(String, Size)>,
    seeds: Vec<(String, f64)>,
    cutoffs: Option<PathBuf>,
    c: Option<usize>,
    runs: Vec<Run>,
    tools: Vec<util::tools::Identity>,
}

impl Preamble {
    /// Take in one `#=` line.
    pub fn absorb(&mut self, row: &toil::MetaRow) -> anyhow::Result<()> {
        match row.key() {
            "query" => self.query = Some(size(row, 0)?),
            "target" => self.targets.push((word(row, 0), size(row, 1)?)),
            "seed" => {
                let wall = row.get(1).context("a `#= seed` line wants a wall time")?;
                self.seeds.push((word(row, 0), wall.parse()?));
            }
            "cutoffs" => {
                let path = row.get(0).context("a `#= cutoffs` line wants a path")?;

                self.cutoffs = Some(PathBuf::from(path));
                self.c = words(row)
                    .skip(1)
                    .find_map(|field| field.strip_prefix("c="))
                    .map(str::parse)
                    .transpose()?;
            }
            "tool" => {
                let Some([Some(name), Some(version), Some(hash)]) = row.exactly::<3>() else {
                    bail!("a `#= tool` line wants a name, a version and a hash");
                };

                self.tools.push(util::tools::Identity {
                    name: name.to_string(),
                    version: version.to_string(),
                    hash: hash.to_string(),
                });
            }
            "run" => self.runs.push(run(row)?),
            other => bail!("unknown `#= {other}` line"),
        }

        Ok(())
    }

    /// The preamble, once the header has been reached.
    pub fn finish(self) -> anyhow::Result<Meta> {
        let meta = Meta {
            query: self.query.context("no `#= query` line")?,
            targets: self.targets,
            seeds: self.seeds,
            cutoffs: self.cutoffs.context("no `#= cutoffs` line")?,
            c: self.c.context("no `c=` on the `#= cutoffs` line")?,
            runs: self.runs,
            tools: self.tools,
        };

        ensure!(!meta.runs.is_empty(), "no `#= run` lines");
        ensure!(!meta.targets.is_empty(), "no `#= target` lines");

        Ok(meta)
    }
}

/// Word `i` of a `#=` line, or empty where it holds the placeholder: a shard
/// with no name, from a pipeline that only ever searched one target.
fn word(row: &toil::MetaRow, i: usize) -> String {
    row.get(i).unwrap_or_default().to_string()
}

/// Every word of a `#=` line after its key.
fn words<'a>(row: &toil::MetaRow<'a>) -> impl Iterator<Item = &'a str> {
    let row = *row;
    (0..row.len()).filter_map(move |i| row.get(i))
}

fn size(row: &toil::MetaRow, from: usize) -> anyhow::Result<Size> {
    let (Some(count), Some(residues), Some(bytes)) =
        (row.get(from), row.get(from + 1), row.get(from + 2))
    else {
        bail!("a size wants a count, a residue count and a byte count");
    };
    ensure!(
        row.len() == from + 3,
        "a size wants a count, a residue count and a byte count"
    );

    Ok(Size {
        count: count.parse()?,
        residues: residues.parse()?,
        bytes: bytes.parse()?,
    })
}

/// One `#= run` line: a name, a tool, a wall time, then whatever settings told
/// this run apart from the others.
fn run(row: &toil::MetaRow) -> anyhow::Result<Run> {
    let (Some(name), Some(tool), Some(wall)) = (row.get(0), row.get(1), row.get(2)) else {
        bail!("a `#= run` line wants a name, a tool and a wall time");
    };

    let mut params: BTreeMap<String, String> = words(row)
        .skip(3)
        .filter_map(|p| p.split_once('='))
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect();

    let seeds = params.remove(util::manifest::SEEDS);

    Ok(Run {
        name: name.to_string(),
        tool: Tool::parse(tool)?,
        wall_s: wall.parse()?,
        params,
        seeds,
    })
}
