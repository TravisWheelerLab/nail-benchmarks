//! Reading back the table a [`pail::Table`] sink wrote.
//!
//! A run records what it did in `manifest.tbl`, and `parse` learns the shape
//! of a pipeline from that rather than from filenames or from knowing which
//! pipeline it is looking at. What makes a row a column of the scores table is
//! a `name` field on it -- see [`Manifest::runs`].
//!
//! The format is whitespace-separated with a `#` header. Two kinds of row never
//! carry fields and so are never runs: a step's own summary line, and the
//! `||` continuation lines a batched step lists its commands on. Both still
//! split cleanly, because the step cell holds a single token either way.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, ensure};

/// The column a run's results file is filed under.
pub const NAME: &str = "name";
/// Which program produced them, which is what says how to read them.
pub const TOOL: &str = "tool";
/// Which target file they were produced against.
pub const SHARD: &str = "shard";
/// What part of the pipeline a command belongs to, for the ones that aren't
/// runs -- seeding, which produces pairs rather than scores.
pub const STAGE: &str = "stage";

/// Everything a pipeline produces lands in one `results/` directory, named
/// `<what>.<shard>`. These three say what the names are.
///
/// They are written by the pipelines and read by `parse`, so they live beside
/// the field names rather than at either end. A pipeline with no shard axis
/// leaves the shard out.
///
/// `seeds` is a reserved stem: a run named it would land on the seed list.
pub fn table_path(results: &Path, name: &str, shard: &str) -> PathBuf {
    results.join(stem(name, shard, Some("tbl")))
}

/// The per-domain table an hmmer run writes alongside its hits.
pub fn dom_path(results: &Path, name: &str, shard: &str) -> PathBuf {
    results.join(stem(name, shard, Some("domtbl")))
}

/// The (query, target) pairs seeding found for one shard. No extension: it is
/// nail's own two-column list, not a hit table, and the bare name says so in a
/// directory of `.tbl`s.
pub fn seeds_path(results: &Path, shard: &str) -> PathBuf {
    results.join(stem("seeds", shard, None))
}

fn stem(name: &str, shard: &str, ext: Option<&str>) -> String {
    let shard = match shard.is_empty() {
        true => String::new(),
        false => format!(".{shard}"),
    };
    let ext = ext.map(|e| format!(".{e}")).unwrap_or_default();

    format!("{name}{shard}{ext}")
}

const WALL: &str = "wall(s)";
const EXIT: &str = "exit";

/// What a batched step puts in the step cell of each of its commands.
const BATCH: &str = "||";

/// Everything one command line of `manifest.tbl` had to say.
pub struct Row {
    cells: BTreeMap<String, String>,
}

impl Row {
    /// One cell, or `None` where the table wrote a dash. Absent and "had
    /// nothing to say" are the same answer here, since a field a command never
    /// set renders as a dash in every row of the block.
    pub fn get(&self, key: &str) -> Option<&str> {
        match self.cells.get(key).map(String::as_str) {
            Some("-") | None => None,
            Some(value) => Some(value),
        }
    }

    pub fn wall_s(&self) -> Option<f64> {
        self.get(WALL)?.parse().ok()
    }

    /// Whether the command finished successfully. A row with no exit code
    /// never ran, which is not success.
    pub fn ok(&self) -> bool {
        self.get(EXIT) == Some("0")
    }

    /// Whether this command ran alongside the others in its step rather than
    /// after them.
    ///
    /// A batched step marks its commands with `||` in the step cell. Their
    /// wall clocks overlap, so the longest of them is what the step took --
    /// adding them up would report the work rather than the time. A step
    /// holding one command collapses to a single row and reads as serial,
    /// which for one command comes to the same thing either way.
    pub fn batched(&self) -> bool {
        self.cells.get("step").map(String::as_str) == Some(BATCH)
    }

    /// Everything this row set that isn't part of the contract: the parameters
    /// that tell one run of a tool apart from another.
    pub fn params(&self) -> BTreeMap<String, String> {
        self.cells
            .iter()
            .filter(|(key, _)| {
                !matches!(key.as_str(), NAME | TOOL | SHARD | STAGE) && !is_metric(key)
            })
            .filter(|(_, value)| value.as_str() != "-")
            .map(|(key, value)| (key.clone(), value.clone()))
            .collect()
    }
}

/// The columns every block ends with, which are the pipeline's own accounting
/// rather than anything a benchmark asked for.
fn is_metric(key: &str) -> bool {
    matches!(
        key,
        "step"
            | "cmd"
            | "wall(s)"
            | "user(s)"
            | "sys(s)"
            | "cpu(%)"
            | "max_rss"
            | "cpus"
            | "exit"
            | "status"
            | "argv"
    )
}

pub struct Manifest {
    rows: Vec<Row>,
}

impl Manifest {
    pub fn read(path: &Path) -> anyhow::Result<Manifest> {
        let text = std::fs::read_to_string(path).with_context(|| {
            format!(
                "failed to read {}; has this pipeline been run?",
                path.display()
            )
        })?;

        let mut names: Vec<String> = Vec::new();
        let mut rows = Vec::new();

        for line in text.lines() {
            if let Some(rest) = line.strip_prefix('#') {
                // the first comment is the header, the second is the rule under it
                let rest = rest.trim_start();
                if names.is_empty() && !rest.starts_with('-') {
                    names = rest.split_whitespace().map(str::to_string).collect();
                }
                continue;
            }

            if line.trim().is_empty() {
                continue;
            }

            // argv is last and full of spaces, so the zip stops at its first
            // word. nothing reads it.
            rows.push(Row {
                cells: names
                    .iter()
                    .cloned()
                    .zip(line.split_whitespace().map(str::to_string))
                    .collect(),
            });
        }

        ensure!(!names.is_empty(), "no header in {}", path.display());
        Ok(Manifest { rows })
    }

    /// The rows that name a run and finished, in the order they were declared.
    ///
    /// A step summary carries no fields and a batch continuation carries the
    /// same fields as the command it belongs to, so naming a run is what marks
    /// one out rather than any guess at the row's shape.
    ///
    /// A command that failed is left out. Its results table is missing or
    /// half-written, and a half-written one is worse than an absent column:
    /// it reads as a run that simply found less.
    pub fn runs(&self) -> impl Iterator<Item = &Row> {
        self.rows
            .iter()
            .filter(|row| row.get(NAME).is_some() && row.ok())
    }

    /// Rows that name a run but did not finish, so a caller can say what is
    /// missing rather than quietly reporting a smaller table.
    pub fn failed(&self) -> impl Iterator<Item = &Row> {
        self.rows
            .iter()
            .filter(|row| row.get(NAME).is_some() && !row.ok())
    }

    /// Command rows belonging to one stage of the pipeline, such as `seed`.
    pub fn stage(&self, stage: &str) -> impl Iterator<Item = &Row> {
        self.rows
            .iter()
            .filter(move |row| row.get(STAGE) == Some(stage))
    }
}
