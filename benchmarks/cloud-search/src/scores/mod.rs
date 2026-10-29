//! One row per query/target pair, and what nail's runs made of it.
//!
//! [`runs`] holds the table, a score per run. What is in this module is what
//! the writer and the reader share: the runs a ledger declares, the query
//! numbering, the cutoffs, and the preamble. [`shard`] reads a shard's results
//! into pairs and [`collect`] turns shards into a file; [`frame`] reads one
//! back.
//!
//! A pair earns a row by clearing some run's cutoff. A family short of any
//! cutoff keeps every pair a run reported, so an analysis that leaves it out
//! can count what it left.

pub mod analyze;
pub mod collect;
pub mod frame;
pub mod parse;
pub mod runs;
mod scan;
pub mod shard;

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};

use libsail::format::Format;
use libsail::index::Reader;
use libsail::seq::p7hmm::leng_of;

use util::ledger::{self, Ledger};
use util::set::Set;

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

/// One column: a named cell of nail at one parameterization, over every rep
/// it ran.
#[derive(Clone, Debug)]
pub struct Run {
    /// The cell's name, without the rep.
    pub name: String,
    /// What each rep cost, summed over every shard it covered, in rep order.
    pub walls: Vec<f64>,
    /// Whatever else the commands recorded: the settings that tell this run
    /// apart from the others.
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
    /// The name rep 1's results tables are written under.
    pub file: String,
    /// In manifest order. Empty string for a pipeline that never named one.
    pub shards: Vec<String>,
}

/// One rep of a cell, as the ledger holds it.
struct Rep {
    rep: usize,
    name: String,
    wall_s: f64,
    shards: Vec<String>,
}

/// The cells a pipeline declared, in the order it declared them, each with
/// every rep it ran.
///
/// # Errors
///
/// Fails naming the cell where a rep reported a different number of hits on
/// some shard than rep 1 did: nail is deterministic, so the reps differ only
/// in what they cost.
pub fn runs(ran: &Ledger, results: &Path) -> anyhow::Result<Vec<Column>> {
    let mut cells: Vec<(Run, Vec<Rep>)> = Vec::new();

    for column in ran.columns()? {
        // a run whose rows are all `*` says what it cost without saying
        // which targets it searched, and a table read by shard has nothing
        // to open. every row of its column would be a dash
        ensure!(
            !column.shards.is_empty(),
            "run {:?} covers no shard of its own in ledger.tbl: every row of it is `{}`",
            column.name,
            ledger::EVERY_SHARD,
        );

        // every table is read as nail's, so another tool's would read
        // as hits in the wrong columns rather than fail
        ensure!(
            column.tool == "nail",
            "run {:?} is {}'s, and this reads nail's tables only",
            column.name,
            column.tool,
        );

        let mut params = column.params;
        let seeds = params.remove(util::manifest::SEEDS);

        // a ledger written before reps has no `rep` setting, and is one rep
        let (cell, rep) = match params.remove(crate::run::REP) {
            Some(rep) => {
                let rep: usize = rep
                    .parse()
                    .with_context(|| format!("run {:?} has rep {rep:?}", column.name))?;
                let cell = column
                    .name
                    .strip_suffix(&format!("-r{rep}"))
                    .with_context(|| {
                        format!(
                            "run {:?} is rep {rep} and its name does not end -r{rep}",
                            column.name
                        )
                    })?;
                (cell.to_string(), rep)
            }
            None => (column.name.clone(), 1),
        };

        let rep = Rep {
            rep,
            name: column.name,
            wall_s: column.wall_s,
            shards: column.shards,
        };

        match cells.iter_mut().find(|(run, _)| run.name == cell) {
            Some((run, reps)) => {
                ensure!(
                    run.params == params,
                    "the reps of {cell} were run with different settings"
                );
                reps.push(rep);
            }
            None => cells.push((
                Run {
                    name: cell,
                    walls: Vec::new(),
                    params,
                    seeds,
                },
                vec![rep],
            )),
        }
    }

    cells
        .into_iter()
        .map(|(mut run, mut reps)| {
            reps.sort_by_key(|rep| rep.rep);

            ensure!(
                reps[0].rep == 1,
                "{} has no rep 1 in the ledger to read its hits from",
                run.name
            );
            for pair in reps.windows(2) {
                ensure!(
                    pair[0].rep != pair[1].rep,
                    "{} has rep {} twice",
                    run.name,
                    pair[0].rep
                );
            }

            let first = &reps[0];
            for rep in &reps[1..] {
                same_hits(results, &run.name, first, rep)?;
            }

            run.walls = reps.iter().map(|rep| rep.wall_s).collect();
            let first = reps.swap_remove(0);

            Ok(Column {
                run,
                file: first.name,
                shards: first.shards,
            })
        })
        .collect()
}

/// Whether `rep` reported as many hits as `first` on every shard.
fn same_hits(results: &Path, cell: &str, first: &Rep, rep: &Rep) -> anyhow::Result<()> {
    let mut want = first.shards.clone();
    let mut got = rep.shards.clone();
    want.sort();
    got.sort();
    ensure!(
        want == got,
        "{cell}: rep {} searched other shards than rep 1",
        rep.rep
    );

    for shard in &first.shards {
        let want = hits(&util::manifest::table_path(results, &first.name, shard))?;
        let got = hits(&util::manifest::table_path(results, &rep.name, shard))?;

        ensure!(
            want == got,
            "{cell}: rep {} reported {got} hits on shard {shard:?} and rep 1 reported {want}",
            rep.rep
        );
    }

    Ok(())
}

/// The rows of a results table: every line but blanks and `#` lines.
fn hits(path: &Path) -> anyhow::Result<usize> {
    use std::io::BufRead;

    let file =
        std::fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;

    let mut n = 0;
    for line in std::io::BufReader::new(file).split(b'\n') {
        let line = line.with_context(|| format!("failed to read {}", path.display()))?;
        n += usize::from(!line.is_empty() && line[0] != b'#');
    }

    Ok(n)
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

/// nail's score threshold on each family, by query id.
pub struct Cutoffs {
    nail: Vec<Option<f32>>,
}

impl Cutoffs {
    /// The threshold a run is held to on one family.
    pub fn get(&self, query: u32) -> Option<f32> {
        self.nail.get(query as usize).copied().flatten()
    }

    /// Whether nail has a threshold on one family.
    pub fn complete(&self, query: u32) -> bool {
        self.get(query).is_some()
    }

    /// One score per family, out of the decoys the calibration scored it
    /// against.
    pub fn read(path: &Path, c: usize, queries: &Queries) -> anyhow::Result<Cutoffs> {
        let mut out = Cutoffs {
            nail: vec![None; queries.len()],
        };

        parse(path, c, |family, nail| {
            // a family the calibration knows and this query set does not is
            // not an error: the calibration covers every family Pfam has
            if let Some(id) = queries.id(family) {
                out.nail[id as usize] = nail;
            }
        })?;

        if out.nail.iter().all(Option::is_none) {
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
pub struct Named(HashMap<String, Option<f32>>);

impl Named {
    /// nail's threshold on one family.
    pub fn get(&self, family: &str) -> Option<f32> {
        self.0.get(family).copied().flatten()
    }

    pub fn read(path: &Path, c: usize) -> anyhow::Result<Named> {
        let mut out = HashMap::new();

        parse(path, c, |family, nail| {
            out.insert(family.to_string(), nail);
        })?;

        ensure!(
            out.values().any(Option::is_some),
            "no usable cutoffs at index {c} in {}",
            path.display()
        );

        Ok(Named(out))
    }
}

/// The families `cutoffs.tbl` gives nail a threshold for at index `c`. Every
/// other family is left out of the analyses.
pub fn complete(path: &Path, c: usize) -> anyhow::Result<HashSet<String>> {
    let mut out = HashSet::new();
    parse(path, c, |family, nail| {
        if nail.is_some() {
            out.insert(family.to_string());
        }
    })?;

    Ok(out)
}

/// Read `cutoffs.tbl`, handing each family nail's threshold at index `c`.
///
/// `cutoffs.tbl` names its columns `<tool>_1..<tool>_5` and `<tool>_n`, so the
/// column to read is found by name rather than by counting -- which is what
/// lets the calibration add a tool without moving anything here.
fn parse(path: &Path, c: usize, mut row: impl FnMut(&str, Option<f32>)) -> anyhow::Result<()> {
    let table =
        toil::Table::read(path).with_context(|| format!("failed to read {}", path.display()))?;

    let name = format!("nail_{}", c + 1);
    let at = table
        .index(&name)
        .with_context(|| format!("{} has no {name} column", path.display()))?;

    for cells in table.rows() {
        let Some(family) = cells.get(0) else {
            continue;
        };

        // a zero is a family nail learned nothing about
        let nail = cells
            .get(at)
            .and_then(|x| x.parse::<f32>().ok())
            .filter(|s| *s > 0.0);

        row(family, nail);
    }

    Ok(())
}

// --------------------------------------------------------------------- meta

/// What a `runs.tbl` says about itself, above the header.
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
    pub fn write<W: std::io::Write>(&self, out: &mut toil::Stream<W>) -> std::io::Result<()> {
        out.meta("format", [runs::FORMAT])?;

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
            let walls: Vec<String> = run.walls.iter().map(|w| format!("{w:.4}")).collect();
            let words = [run.name.clone(), walls.join(",")]
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

/// One `#= run` line: a name, a wall time per rep joined by commas, then
/// whatever settings told this run apart from the others.
fn run(row: &toil::MetaRow) -> anyhow::Result<Run> {
    let (Some(name), Some(wall)) = (row.get(0), row.get(1)) else {
        bail!("a `#= run` line wants a name and a wall time");
    };

    let mut params: BTreeMap<String, String> = words(row)
        .skip(2)
        .filter_map(|p| p.split_once('='))
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect();

    let seeds = params.remove(util::manifest::SEEDS);

    Ok(Run {
        name: name.to_string(),
        walls: wall
            .split(',')
            .map(str::parse)
            .collect::<Result<_, _>>()
            .with_context(|| format!("run {name:?} has wall times {wall:?}"))?,
        params,
        seeds,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A ledger of one cell run twice on shard `s`, and the results directory
    /// holding its tables, with rep 2's carrying `second` hits to rep 1's two.
    fn twice(name: &str, second: usize) -> (Ledger, PathBuf) {
        let results =
            std::env::temp_dir().join(format!("cloud-search-reps-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&results).unwrap();

        let table = |hits: usize| {
            let mut text = String::from("# query target score\n");
            for i in 0..hits {
                text.push_str(&format!("alpha MGYP{i:012} 30.0\n"));
            }
            text
        };
        let r1 = util::manifest::table_path(&results, "A2.0-B4.0-r1", "s");
        let r2 = util::manifest::table_path(&results, "A2.0-B4.0-r2", "s");
        std::fs::write(r1, table(2)).unwrap();
        std::fs::write(r2, table(second)).unwrap();

        let row = |rep: usize, wall: f64| ledger::Row {
            name: format!("A2.0-B4.0-r{rep}"),
            tool: "nail".to_string(),
            shard: "s".to_string(),
            stage: String::new(),
            params: BTreeMap::from([
                ("A".to_string(), "2.0".to_string()),
                ("B".to_string(), "4.0".to_string()),
                (crate::run::REP.to_string(), rep.to_string()),
            ]),
            wall_s: Some(wall),
            cpu_s: None,
            max_rss_kb: None,
        };

        (Ledger::from_rows(vec![row(2, 5.0), row(1, 3.0)]), results)
    }

    #[test]
    fn reps_fold_into_one_column_read_from_rep_1() {
        let (ran, results) = twice("same", 2);
        let columns = runs(&ran, &results).unwrap();

        assert_eq!(columns.len(), 1);
        assert_eq!(columns[0].run.name, "A2.0-B4.0");
        assert_eq!(columns[0].file, "A2.0-B4.0-r1");
        assert_eq!(columns[0].run.walls, [3.0, 5.0]);
        assert!(!columns[0].run.params.contains_key(crate::run::REP));
    }

    #[test]
    fn a_rep_with_other_hits_than_rep_1_fails_naming_the_cell() {
        let (ran, results) = twice("differ", 3);
        let Err(err) = runs(&ran, &results) else {
            panic!("rep 2 found a hit rep 1 did not, and nothing failed");
        };
        let err = err.to_string();

        assert!(err.contains("A2.0-B4.0: rep 2 reported 3 hits"), "{err}");
    }
}
