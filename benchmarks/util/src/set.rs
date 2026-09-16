//! What a build produced, described so a search does not have to know who
//! built it.
//!
//! A pipeline needs four things of an input set: what to search with, what to
//! search, which representations of the query exist, and whatever tells one
//! unit of work from another. `set.tbl` holds exactly those, one row per
//! search unit, and the builder writes it in the same pass that writes the
//! files.
//!
//! That is the same bargain [`crate::ledger`] strikes one seam later. A ledger
//! row is a spine plus an open map of settings, so an analysis reads a run's
//! shape out of the table rather than out of the filenames; a set row is a
//! spine plus an open map of attributes, so a search reads a set's shape the
//! same way. A recipe that deals shards and a recipe that nests rungs produce
//! the same table with different attribute columns, and a pipeline reading it
//! does not learn which ran.
//!
//! ```text
//! # unit query_hmm       query_sto       query_db          target        shard residues
//! # ---- --------------- --------------- ----------------- ------------- ----- --------
//!   1    queries/q.hmm   queries/q.sto   queries/qDB/qDB   targets/1.fa  1     34129933
//! ```
//!
//! Paths are relative to the set's own directory, so a set is one tree that
//! can be moved or linked without rewriting the table. A cell left empty means
//! the builder did not produce that representation, and asking for it fails
//! naming the set rather than failing later on a path that was never there.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};

use crate::tbl;

/// What the manifest is called, in a set's directory.
pub const FILE: &str = "set.tbl";

const UNIT: &str = "unit";
const QUERY_HMM: &str = "query_hmm";
const QUERY_STO: &str = "query_sto";
const QUERY_FA: &str = "query_fa";
const QUERY_DB: &str = "query_db";
const TARGET: &str = "target";

/// The columns every set has, in the order they are written.
const SPINE: [&str; 6] = [UNIT, QUERY_HMM, QUERY_STO, QUERY_FA, QUERY_DB, TARGET];

/// The meta key a builder stamps its shape under.
const SHAPE: &str = "shape";

// ---

/// One of the columns a unit can carry a path in.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Rep {
    QueryHmm,
    QuerySto,
    QueryFa,
    QueryDb,
    Target,
}

impl Rep {
    /// The column this is, which is also what an error calls it.
    pub fn column(self) -> &'static str {
        match self {
            Rep::QueryHmm => QUERY_HMM,
            Rep::QuerySto => QUERY_STO,
            Rep::QueryFa => QUERY_FA,
            Rep::QueryDb => QUERY_DB,
            Rep::Target => TARGET,
        }
    }
}

/// What a benchmark needs of a set, as a thing that can be checked rather than
/// discovered one missing path at a time.
///
/// The open half of `set.tbl` is what lets one manifest describe a deal of
/// shards and a nest of rungs. The cost of that is that a set and the benchmark
/// reading it can disagree without anything saying so: a ladder loads fine
/// against a sweep expecting shards, and the sweep searches the product of both
/// axes as if it were a list of targets. A shape is the declaration that closes
/// that -- the builder stamps one, the benchmark names the one it wants, and
/// [`Set::check`] holds them to each other before a search starts.
pub struct Shape {
    pub name: &'static str,
    /// The columns every unit must carry a path in.
    pub needs: &'static [Rep],
    /// The attributes every unit must carry, whatever their values.
    pub attrs: &'static [&'static str],
}

/// The shapes the benchmarks in this repo build and read.
///
/// They live here rather than with the benchmarks because a shape is the
/// agreement between a builder and a reader, and neither side owns it.
pub mod shape {
    use super::{Rep, Shape};

    /// One query set against target shards of equal size. The shard is a unit
    /// of work rather than a variable, so every unit names the same query.
    pub const FIXED: Shape = Shape {
        name: "fixed",
        needs: &[Rep::QueryHmm, Rep::QuerySto, Rep::QueryDb, Rep::Target],
        // seqs, residues and bytes are read by the analyses rather than by a
        // search: a shape is what a whole benchmark needs, not what its first
        // stage opens
        attrs: &["shard", "seqs", "residues", "bytes"],
    };

    /// Nested rungs on both axes, each a prefix of the one above.
    pub const LADDER: Shape = Shape {
        name: "ladder",
        needs: &[Rep::QueryHmm, Rep::QuerySto, Rep::QueryDb, Rep::Target],
        attrs: &[
            "query_rung",
            "target_rung",
            "query_residues",
            "target_residues",
        ],
    };

    /// One query against one target, paired rather than crossed.
    pub const PAIRS: Shape = Shape {
        name: "pairs",
        needs: &[Rep::QueryFa, Rep::Target],
        attrs: &["pair", "query_residues", "residues"],
    };

    /// Per-family decoys, forward and reversed, drawn from another set.
    pub const DECOYS: Shape = Shape {
        name: "decoys",
        needs: &[Rep::QueryHmm, Rep::QuerySto, Rep::Target],
        attrs: &["family", "direction"],
    };
}

// ---

/// One search a tool will be asked to run: a query against a target.
pub struct Row {
    /// What tells this unit from the others, and what a result table is named
    /// after.
    pub unit: String,
    /// The query as profiles, relative to the set root, or empty.
    pub query_hmm: String,
    /// The query as alignments, relative to the set root, or empty.
    pub query_sto: String,
    /// The query as plain sequences, relative to the set root, or empty. What
    /// a tool searched in sequence mode reads, where a profile search reads
    /// `query_hmm`.
    pub query_fa: String,
    /// The query as an mmseqs profile database, relative to the set root, or
    /// empty.
    pub query_db: String,
    /// What is being searched, relative to the set root.
    pub target: String,
    /// Everything else a recipe wrote down: a shard number, a rung, a family,
    /// a residue count.
    pub attrs: BTreeMap<String, String>,
}

impl Row {
    /// A row with only the two things every set has.
    pub fn new(unit: impl Into<String>, target: impl Into<String>) -> Row {
        Row {
            unit: unit.into(),
            query_hmm: String::new(),
            query_sto: String::new(),
            query_fa: String::new(),
            query_db: String::new(),
            target: target.into(),
            attrs: BTreeMap::new(),
        }
    }

    pub fn query_hmm(mut self, path: impl Into<String>) -> Row {
        self.query_hmm = path.into();
        self
    }

    pub fn query_sto(mut self, path: impl Into<String>) -> Row {
        self.query_sto = path.into();
        self
    }

    pub fn query_fa(mut self, path: impl Into<String>) -> Row {
        self.query_fa = path.into();
        self
    }

    pub fn query_db(mut self, path: impl Into<String>) -> Row {
        self.query_db = path.into();
        self
    }

    pub fn attr(mut self, key: impl Into<String>, value: impl std::fmt::Display) -> Row {
        self.attrs.insert(key.into(), value.to_string());
        self
    }

    fn cell(&self, rep: Rep) -> &str {
        match rep {
            Rep::QueryHmm => &self.query_hmm,
            Rep::QuerySto => &self.query_sto,
            Rep::QueryFa => &self.query_fa,
            Rep::QueryDb => &self.query_db,
            Rep::Target => &self.target,
        }
    }
}

// ---

/// A set on disk: its directory, its rows, and what it says about itself.
pub struct Set {
    root: PathBuf,
    rows: Vec<Row>,
    /// Provenance, as written to and read from the table's `#=` lines.
    meta: BTreeMap<String, String>,
}

impl Set {
    pub fn new(root: impl Into<PathBuf>, rows: Vec<Row>) -> Set {
        Set {
            root: root.into(),
            rows,
            meta: BTreeMap::new(),
        }
    }

    /// Record where this set came from: the recipe that ran, its parameters,
    /// the sources it drew from, or the set it derives from.
    pub fn says(mut self, key: impl Into<String>, value: impl std::fmt::Display) -> Set {
        self.meta.insert(key.into(), value.to_string());
        self
    }

    /// The set in a directory, read off its own `set.tbl`.
    pub fn load(dir: impl AsRef<Path>) -> anyhow::Result<Set> {
        let dir = dir.as_ref();
        let path = dir.join(FILE);

        ensure!(
            path.is_file(),
            "no {} in {}; has the set been built?",
            FILE,
            dir.display()
        );

        Set::read(dir, &path)
    }

    /// The set in a directory, checked against the shape the caller needs.
    ///
    /// The check is here rather than at the first path a search opens, so a set
    /// of the wrong shape fails before any tool runs.
    pub fn load_as(dir: impl AsRef<Path>, shape: &Shape) -> anyhow::Result<Set> {
        let set = Set::load(dir)?;
        set.check(shape)?;
        Ok(set)
    }

    /// Hold this set to a shape: what it says it is, then what every unit
    /// carries.
    pub fn check(&self, shape: &Shape) -> anyhow::Result<()> {
        let at = self.root.display();

        if let Some(said) = self.said(SHAPE) {
            ensure!(
                said == shape.name,
                "the set at {at} was built as {said:?}, but this reads a {:?} set",
                shape.name
            );
        }

        for (row, unit) in self.rows.iter().zip(self.units()) {
            for &rep in shape.needs {
                ensure!(
                    !row.cell(rep).is_empty(),
                    "the set at {at} is not a {:?} set: unit {:?} has no {}",
                    shape.name,
                    unit.name(),
                    rep.column()
                );
            }

            for attr in shape.attrs {
                ensure!(
                    row.attrs.contains_key(*attr),
                    "the set at {at} is not a {:?} set: unit {:?} has no {attr}",
                    shape.name,
                    unit.name()
                );
            }
        }

        Ok(())
    }

    pub fn read(root: impl Into<PathBuf>, path: &Path) -> anyhow::Result<Set> {
        let root = root.into();
        let table = tbl::read(path)?;

        let attrs: Vec<&String> = table
            .headers
            .iter()
            .filter(|h| !SPINE.contains(&h.as_str()))
            .collect();

        let mut rows = Vec::new();
        for cells in &table.cells {
            let cell = |key: &str| match cells.get(key).map(String::as_str) {
                Some("-") | None => "",
                Some(value) => value,
            };

            let target = cell(TARGET);
            ensure!(
                !target.is_empty(),
                "{} has a row with no {TARGET}",
                path.display()
            );

            rows.push(Row {
                unit: cell(UNIT).to_string(),
                query_hmm: cell(QUERY_HMM).to_string(),
                query_sto: cell(QUERY_STO).to_string(),
                query_fa: cell(QUERY_FA).to_string(),
                query_db: cell(QUERY_DB).to_string(),
                target: target.to_string(),
                attrs: attrs
                    .iter()
                    .filter(|key| !cell(key).is_empty())
                    .map(|key| ((*key).clone(), cell(key).to_string()))
                    .collect(),
            });
        }

        ensure!(!rows.is_empty(), "{} names no units", path.display());

        Ok(Set {
            root,
            rows,
            meta: read_meta(&table.meta),
        })
    }

    /// Write the table into the set's own directory, which is where every
    /// reader looks for it.
    pub fn save(&self) -> anyhow::Result<()> {
        self.write(&self.root.join(FILE))
    }

    pub fn write(&self, path: &Path) -> anyhow::Result<()> {
        let keys: Vec<String> = {
            let mut keys: Vec<String> = self
                .rows
                .iter()
                .flat_map(|row| row.attrs.keys().cloned())
                .collect();
            keys.sort_unstable();
            keys.dedup();
            keys
        };

        let mut headers: Vec<String> = SPINE.iter().map(|h| h.to_string()).collect();
        headers.extend(keys.iter().cloned());

        let rows: Vec<Vec<String>> = self
            .rows
            .iter()
            .map(|row| {
                let mut cells = vec![
                    dash(&row.unit),
                    dash(&row.query_hmm),
                    dash(&row.query_sto),
                    dash(&row.query_fa),
                    dash(&row.query_db),
                    dash(&row.target),
                ];
                cells.extend(keys.iter().map(|key| match row.attrs.get(key) {
                    Some(value) => value.clone(),
                    None => "-".to_string(),
                }));
                cells
            })
            .collect();

        let meta: String = self
            .meta
            .iter()
            .map(|(key, value)| format!("#= {key} {value}\n"))
            .collect();

        tbl::write(
            path,
            tbl::Table {
                meta: &meta,
                headers: &headers,
                rows: &rows,
                ragged_last: false,
            },
        )
    }

    pub fn root(&self) -> &Path {
        &self.root
    }

    pub fn len(&self) -> usize {
        self.rows.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// What the set says about where it came from.
    pub fn said(&self, key: &str) -> Option<&str> {
        self.meta.get(key).map(String::as_str)
    }

    /// Every unit, with its paths already resolved against the set root.
    pub fn units(&self) -> impl Iterator<Item = Unit<'_>> {
        self.rows.iter().map(|row| Unit {
            row,
            root: &self.root,
        })
    }

    /// The units whose `key` attribute is `value`, for a pipeline that walks
    /// one axis at a time.
    pub fn where_attr<'a>(&'a self, key: &'a str, value: &'a str) -> impl Iterator<Item = Unit<'a>> {
        self.units().filter(move |u| u.attr(key) == Some(value))
    }

    /// Every distinct value of `key`, in the order the rows give them.
    pub fn values(&self, key: &str) -> Vec<String> {
        let mut seen = Vec::new();
        for row in &self.rows {
            if let Some(value) = row.attrs.get(key) {
                if !seen.iter().any(|v| v == value) {
                    seen.push(value.clone());
                }
            }
        }
        seen
    }
}

// ---

/// One row with the set root it is relative to, so a caller gets paths rather
/// than strings.
#[derive(Clone, Copy)]
pub struct Unit<'a> {
    row: &'a Row,
    root: &'a Path,
}

impl<'a> Unit<'a> {
    pub fn name(&self) -> &'a str {
        &self.row.unit
    }

    pub fn query_hmm(&self) -> anyhow::Result<PathBuf> {
        self.resolve(QUERY_HMM, &self.row.query_hmm)
    }

    pub fn query_sto(&self) -> anyhow::Result<PathBuf> {
        self.resolve(QUERY_STO, &self.row.query_sto)
    }

    pub fn query_fa(&self) -> anyhow::Result<PathBuf> {
        self.resolve(QUERY_FA, &self.row.query_fa)
    }

    pub fn query_db(&self) -> anyhow::Result<PathBuf> {
        self.resolve(QUERY_DB, &self.row.query_db)
    }

    pub fn target(&self) -> anyhow::Result<PathBuf> {
        self.resolve(TARGET, &self.row.target)
    }

    pub fn attr(&self, key: &str) -> Option<&'a str> {
        self.row.attrs.get(key).map(String::as_str)
    }

    /// An attribute the caller cannot do without, named in the error rather
    /// than unwrapped.
    pub fn need(&self, key: &str) -> anyhow::Result<&'a str> {
        self.attr(key)
            .with_context(|| format!("unit {:?} has no {key}", self.row.unit))
    }

    pub fn number(&self, key: &str) -> anyhow::Result<u64> {
        let text = self.need(key)?;
        text.parse()
            .with_context(|| format!("unit {:?} has a {key} that is not a number: {text:?}", self.row.unit))
    }

    fn resolve(&self, what: &str, cell: &str) -> anyhow::Result<PathBuf> {
        if cell.is_empty() {
            bail!(
                "the set at {} has no {what}, which unit {:?} needs",
                self.root.display(),
                self.row.unit
            );
        }

        Ok(self.root.join(cell))
    }
}

// ---

fn dash(cell: &str) -> String {
    match cell.is_empty() {
        true => "-".to_string(),
        false => cell.to_string(),
    }
}

fn read_meta(lines: &[String]) -> BTreeMap<String, String> {
    lines
        .iter()
        .filter_map(|line| line.strip_prefix("#="))
        .filter_map(|rest| {
            let mut parts = rest.trim().splitn(2, char::is_whitespace);
            let key = parts.next()?;
            Some((key.to_string(), parts.next().unwrap_or("").trim().to_string()))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set(root: &str) -> Set {
        Set::new(
            root,
            vec![
                Row::new("1", "targets/1.fa")
                    .query_hmm("queries/query.hmm")
                    .attr("shard", 1)
                    .attr("residues", 34129933u64),
                Row::new("2", "targets/2.fa")
                    .query_hmm("queries/query.hmm")
                    .attr("shard", 2)
                    .attr("residues", 34_000_000u64),
            ],
        )
        .says("recipe", "fixed")
    }

    // one directory per test: these run in parallel in one process, and they
    // all write a file called set.tbl
    fn round_trip(set: &Set, what: &str) -> Set {
        let dir =
            std::env::temp_dir().join(format!("util-set-{}-{what}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let path = dir.join(FILE);
        set.write(&path).unwrap();
        let back = Set::read(set.root(), &path).unwrap();

        std::fs::remove_dir_all(&dir).ok();
        back
    }

    #[test]
    fn a_set_survives_a_round_trip() {
        let before = set("/tmp/x");
        let after = round_trip(&before, "trip");

        assert_eq!(after.len(), 2);
        assert_eq!(after.said("recipe"), Some("fixed"));

        let units: Vec<Unit<'_>> = after.units().collect();
        assert_eq!(units[0].name(), "1");
        assert_eq!(units[0].target().unwrap(), PathBuf::from("/tmp/x/targets/1.fa"));
        assert_eq!(
            units[0].query_hmm().unwrap(),
            PathBuf::from("/tmp/x/queries/query.hmm")
        );
    }

    #[test]
    fn unknown_columns_survive_a_round_trip() {
        let after = round_trip(&set("/tmp/x"), "columns");
        let units: Vec<Unit<'_>> = after.units().collect();

        assert_eq!(units[0].attr("shard"), Some("1"));
        assert_eq!(units[0].number("residues").unwrap(), 34129933);
        assert_eq!(units[1].attr("shard"), Some("2"));
    }

    #[test]
    fn a_representation_that_was_not_built_names_the_set() {
        let after = round_trip(&set("/tmp/x"), "absent");
        let unit = after.units().next().unwrap();

        let err = unit.query_db().unwrap_err().to_string();
        assert!(err.contains("query_db"), "{err}");
        assert!(err.contains("/tmp/x"), "{err}");
    }

    #[test]
    fn the_values_of_an_attribute_come_back_in_row_order() {
        let after = round_trip(&set("/tmp/x"), "values");
        assert_eq!(after.values("shard"), vec!["1", "2"]);
    }

    /// A set that satisfies `shape::FIXED`, for the checks to work against.
    fn fixed_set() -> Set {
        Set::new(
            "/tmp/x",
            (1..=2)
                .map(|shard| {
                    Row::new(shard.to_string(), format!("targets/{shard}.fa"))
                        .query_hmm("queries/query.hmm")
                        .query_sto("queries/query.sto")
                        .query_db("queries/queryDB/queryDB")
                        .attr("shard", shard)
                        .attr("seqs", 200)
                        .attr("residues", 41677u64)
                        .attr("bytes", 47340u64)
                })
                .collect(),
        )
        .says("shape", "fixed")
    }

    #[test]
    fn a_set_that_satisfies_a_shape_passes() {
        assert!(fixed_set().check(&shape::FIXED).is_ok());
    }

    #[test]
    fn a_missing_representation_names_the_unit_and_the_column() {
        let mut set = fixed_set();
        set.rows[1].query_db.clear();

        let err = set.check(&shape::FIXED).unwrap_err().to_string();
        assert!(err.contains("query_db"), "{err}");
        assert!(err.contains("\"2\""), "{err}");
    }

    #[test]
    fn a_missing_attribute_names_the_unit_and_the_attribute() {
        let mut set = fixed_set();
        set.rows[0].attrs.remove("residues");

        let err = set.check(&shape::FIXED).unwrap_err().to_string();
        assert!(err.contains("residues"), "{err}");
        assert!(err.contains("\"1\""), "{err}");
    }

    /// The failure this whole idea exists for: a ladder loads against a sweep
    /// expecting shards, and without the check the sweep searches the product
    /// of both axes as if it were a list of targets.
    #[test]
    fn a_ladder_is_refused_where_a_fixed_set_is_wanted() {
        let ladder = Set::new(
            "/tmp/x",
            vec![
                Row::new("q10.t1000", "targets/1000.fa")
                    .query_hmm("queries/10/query.hmm")
                    .query_sto("queries/10/query.sto")
                    .query_db("queries/10/queryDB/queryDB")
                    .attr("query_rung", 10)
                    .attr("target_rung", 1000)
                    .attr("query_residues", 1175u64)
                    .attr("target_residues", 41677u64),
            ],
        )
        .says("shape", "ladder");

        assert!(ladder.check(&shape::LADDER).is_ok());

        let err = ladder.check(&shape::FIXED).unwrap_err().to_string();
        assert!(err.contains("ladder"), "{err}");
        assert!(err.contains("fixed"), "{err}");
    }

    /// A set built before shapes existed says nothing about itself, so the
    /// columns are all there is to go on.
    #[test]
    fn a_set_that_declares_no_shape_is_checked_on_its_columns() {
        let mut set = fixed_set();
        set.meta.remove("shape");

        assert!(set.check(&shape::FIXED).is_ok());
        assert!(set.check(&shape::PAIRS).is_err());
    }

    #[test]
    fn a_missing_attribute_is_named_rather_than_unwrapped() {
        let after = round_trip(&set("/tmp/x"), "need");
        let err = after.units().next().unwrap().need("rung").unwrap_err().to_string();
        assert!(err.contains("rung"), "{err}");
    }
}
