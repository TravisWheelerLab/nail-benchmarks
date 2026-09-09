//! What a pipeline ran, and what each search cost.
//!
//! An analysis needs four things of a run: what to call it, how to read its
//! table, which target it searched, and how long it took. `ledger.tbl` holds
//! exactly those, one row per run per shard, and two things can write it:
//! [`Ledger::distill`] reads them out of the `manifest.tbl` a `michi` pipeline
//! left behind, and a benchmark that took its results from somewhere else
//! reads them out of whatever `.time` files came back with them.
//!
//! That is the whole reason this sits between the two. A manifest records what
//! only the harness that ran the commands can know -- an exit code, an argv, a
//! step's batching -- so a directory of results from a cluster cannot be
//! written as one without inventing those fields. It can be written as this.
//!
//! The folding happens here rather than in an analysis: a run is one row per
//! shard with its seconds already totalled, so nothing downstream needs to
//! know what a batched step is.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, ensure};

use crate::manifest::{self, Manifest};
use crate::tbl;

/// The shard cell of a row whose wall clock covers every shard of its run at
/// once, for a search timed as a whole somewhere else.
pub const EVERY_SHARD: &str = "*";

const WALL: &str = "wall(s)";

/// What the ledger is called, in a pipeline directory.
pub const FILE: &str = "ledger.tbl";

// ---

/// One run against one shard, or one stage of the pipeline against one shard.
pub struct Row {
    /// The run this belongs to, or empty on a stage row.
    pub name: String,
    /// Which program produced the table, which is what says how to read it.
    /// Empty on a stage row.
    pub tool: String,
    /// Which target it searched, empty for a pipeline with no shard axis, or
    /// [`EVERY_SHARD`].
    pub shard: String,
    /// The part of the pipeline this is, for the commands that are not runs.
    /// Empty on a run row.
    pub stage: String,
    /// What tells one run of a tool apart from another.
    pub params: BTreeMap<String, String>,
    /// Seconds, or `None` for a row that says a run covered this shard
    /// without saying what that cost.
    pub wall_s: Option<f64>,
}

impl Row {
    fn is_run(&self) -> bool {
        self.stage.is_empty()
    }
}

/// One run folded across the shards it covered: what every analysis here
/// wants a column to be.
pub struct Column {
    pub name: String,
    pub tool: String,
    pub params: BTreeMap<String, String>,
    /// Seconds over every shard, or the whole-run figure where one shard's
    /// row covers them all.
    pub wall_s: f64,
    /// The shards it covered, in the order the pipeline named them.
    pub shards: Vec<String>,
}

pub struct Ledger {
    rows: Vec<Row>,
    /// What was left out, as (run or stage, shard), so a caller can say what
    /// is missing rather than quietly reporting a smaller table.
    failed: Vec<(String, String)>,
}

impl Ledger {
    /// The ledger in a pipeline directory, or the manifest distilled where
    /// there is no ledger.
    //
    // a pipeline here writes its ledger when it finishes, and clears the
    // previous one before it starts, so a directory holding a manifest and no
    // ledger is a run that died partway: distilling is what reports it. a
    // directory with no manifest at all is results from somewhere else, where
    // the ledger is the only record there is
    pub fn load(dir: &Path) -> anyhow::Result<Ledger> {
        let ledger = path(dir);
        if ledger.is_file() {
            return Ledger::read(&ledger);
        }

        let manifest = dir.join("manifest.tbl");
        ensure!(
            manifest.is_file(),
            "no {} and no manifest.tbl in {}; has this been run?",
            FILE,
            dir.display()
        );

        Ledger::distill(&manifest)
    }

    pub fn read(path: &Path) -> anyhow::Result<Ledger> {
        let table = tbl::read(path)?;
        let params: Vec<&String> = table
            .headers
            .iter()
            .filter(|h| !is_column(h))
            .collect::<Vec<_>>();

        let mut rows = Vec::new();
        for cells in &table.cells {
            let cell = |key: &str| match cells.get(key).map(String::as_str) {
                Some("-") | None => "",
                Some(value) => value,
            };

            let wall = cell(WALL);
            rows.push(Row {
                name: cell(manifest::NAME).to_string(),
                tool: cell(manifest::TOOL).to_string(),
                shard: cell(manifest::SHARD).to_string(),
                stage: cell(manifest::STAGE).to_string(),
                params: params
                    .iter()
                    .filter(|key| !cell(key).is_empty())
                    .map(|key| ((*key).clone(), cell(key).to_string()))
                    .collect(),
                wall_s: match wall.is_empty() {
                    true => None,
                    false => Some(
                        wall.parse()
                            .with_context(|| format!("{wall:?} is not a wall clock time"))?,
                    ),
                },
            });
        }

        Ok(Ledger {
            rows,
            failed: failed_meta(&table.meta),
        })
    }

    pub fn write(&self, path: &Path) -> anyhow::Result<()> {
        let keys: Vec<String> = {
            let mut keys: Vec<String> = self
                .rows
                .iter()
                .flat_map(|row| row.params.keys().cloned())
                .collect();
            keys.sort_unstable();
            keys.dedup();
            keys
        };

        let mut headers = vec![
            manifest::NAME.to_string(),
            manifest::TOOL.to_string(),
            manifest::SHARD.to_string(),
            manifest::STAGE.to_string(),
        ];
        headers.extend(keys.iter().cloned());
        headers.push(WALL.to_string());

        let rows: Vec<Vec<String>> = self
            .rows
            .iter()
            .map(|row| {
                let mut cells = vec![
                    dash(&row.name),
                    dash(&row.tool),
                    dash(&row.shard),
                    dash(&row.stage),
                ];
                cells.extend(keys.iter().map(|key| match row.params.get(key) {
                    Some(value) => value.clone(),
                    None => "-".to_string(),
                }));
                cells.push(match row.wall_s {
                    Some(wall) => format!("{wall:.2}"),
                    None => "-".to_string(),
                });
                cells
            })
            .collect();

        let meta: String = self
            .failed
            .iter()
            .map(|(what, shard)| format!("#= failed {} {}\n", dash(what), dash(shard)))
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

    pub fn runs(&self) -> impl Iterator<Item = &Row> {
        self.rows.iter().filter(|row| row.is_run())
    }

    /// The rows of one stage of the pipeline, such as `seed`.
    pub fn stage(&self, stage: &str) -> impl Iterator<Item = &Row> {
        self.rows.iter().filter(move |row| row.stage == stage)
    }

    /// The runs, in the order the pipeline declared them.
    ///
    /// A name that turns up against several shards is one column covering all
    /// of them: what changed between those searches is the target rather than
    /// the parameterization.
    pub fn columns(&self) -> anyhow::Result<Vec<Column>> {
        let mut out: Vec<Column> = Vec::new();
        let mut at: BTreeMap<&str, usize> = BTreeMap::new();
        // a whole-run figure replaces the sum rather than joining it
        let mut every: Vec<Option<f64>> = Vec::new();

        for row in self.runs() {
            let i = match at.get(row.name.as_str()) {
                Some(&i) => i,
                None => {
                    ensure!(!row.tool.is_empty(), "run {:?} has no tool", row.name);

                    at.insert(&row.name, out.len());
                    every.push(None);
                    out.push(Column {
                        name: row.name.clone(),
                        tool: row.tool.clone(),
                        params: row.params.clone(),
                        wall_s: 0.0,
                        shards: Vec::new(),
                    });
                    out.len() - 1
                }
            };

            match row.shard == EVERY_SHARD {
                true => every[i] = row.wall_s,
                false => {
                    if !out[i].shards.contains(&row.shard) {
                        out[i].shards.push(row.shard.clone());
                    }
                    out[i].wall_s += row.wall_s.unwrap_or(0.0);
                }
            }
        }

        for (column, every) in out.iter_mut().zip(&every) {
            if let Some(wall) = every {
                column.wall_s = *wall;
            }
        }

        Ok(out)
    }

    /// Every shard any run covered, in the order the runs named them.
    pub fn shards(&self) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();

        for row in self.runs() {
            if row.shard != EVERY_SHARD && !out.contains(&row.shard) {
                out.push(row.shard.clone());
            }
        }

        out
    }

    pub fn failed(&self) -> &[(String, String)] {
        &self.failed
    }

    pub fn from_rows(rows: Vec<Row>) -> Ledger {
        Ledger {
            rows,
            failed: Vec::new(),
        }
    }
}

fn is_column(header: &str) -> bool {
    matches!(
        header,
        manifest::NAME | manifest::TOOL | manifest::SHARD | manifest::STAGE | WALL
    )
}

fn dash(cell: &str) -> String {
    match cell.is_empty() {
        true => "-".to_string(),
        false => cell.to_string(),
    }
}

fn failed_meta(meta: &[String]) -> Vec<(String, String)> {
    let mut out = Vec::new();

    for line in meta {
        let mut fields = line.split_whitespace();
        if fields.next() != Some("#=") || fields.next() != Some("failed") {
            continue;
        }

        let what = fields.next().unwrap_or("-");
        let shard = fields.next().unwrap_or("-");
        out.push((undash(what), undash(shard)));
    }

    out
}

fn undash(cell: &str) -> String {
    match cell {
        "-" => String::new(),
        cell => cell.to_string(),
    }
}

// ------------------------------------------------------------ from a manifest

impl Ledger {
    /// What a `michi` pipeline recorded, folded to one row per run per shard.
    pub fn distill(path: &Path) -> anyhow::Result<Ledger> {
        let manifest = Manifest::read(path)?;
        let mut groups: Vec<Group> = Vec::new();
        let mut at: BTreeMap<(String, String), usize> = BTreeMap::new();

        // michi puts a label in the step cell on a step's own summary row and
        // `|` or `||` on the commands under it, so a cell that is neither
        // starts a step. the commands of one batched step overlap and the
        // longest is what the step took; steps themselves run in sequence
        let mut step = 0;

        for row in manifest.rows() {
            match row.step() {
                Some("|") | Some("||") => {}
                _ => step += 1,
            }

            let (key, stage) = match (row.get(manifest::NAME), row.get(manifest::STAGE)) {
                (Some(name), _) => (name.to_string(), false),
                (None, Some(stage)) => (stage.to_string(), true),
                // a step summary, or the scaffolding around a search: neither
                // names a run nor belongs to a stage
                (None, None) => continue,
            };
            let shard = row.get(manifest::SHARD).unwrap_or_default().to_string();

            let i = *at.entry((key.clone(), shard.clone())).or_insert_with(|| {
                groups.push(Group::new(&key, &shard, stage));
                groups.len() - 1
            });
            groups[i].add(step, row)?;
        }

        let failed: Vec<(String, String)> = groups
            .iter()
            .filter(|group| group.failed)
            .map(|group| (group.what.clone(), group.shard.clone()))
            .collect();

        let mut rows: Vec<Row> = Vec::new();
        for group in groups.iter().filter(|group| !group.stage && group.ran) {
            rows.push(group.row());
        }
        for group in groups.iter().filter(|group| group.stage && group.ran) {
            rows.push(group.row());
        }

        Ok(Ledger { rows, failed })
    }
}

/// The commands of one run, or one stage, against one shard.
struct Group {
    what: String,
    shard: String,
    stage: bool,
    tool: String,
    params: BTreeMap<String, String>,
    /// Seconds per step: the serial commands added up, and the longest of the
    /// batched ones.
    steps: BTreeMap<usize, (f64, f64)>,
    /// Whether any command of this group finished. A group of nothing but
    /// failures is not a run.
    ran: bool,
    failed: bool,
}

impl Group {
    fn new(what: &str, shard: &str, stage: bool) -> Group {
        Group {
            what: what.to_string(),
            shard: shard.to_string(),
            stage,
            tool: String::new(),
            params: BTreeMap::new(),
            steps: BTreeMap::new(),
            ran: false,
            failed: false,
        }
    }

    fn add(&mut self, step: usize, row: &manifest::Row) -> anyhow::Result<()> {
        // a command that failed is left out: its table is missing or
        // half-written, and a half-written one reads as a run that simply
        // found less
        if !row.ok() {
            self.failed = true;
            return Ok(());
        }
        self.ran = true;

        if let Some(tool) = row.get(manifest::TOOL) {
            self.tool = tool.to_string();
        }

        for (key, value) in row.params() {
            if let Some(had) = self.params.get(&key) {
                ensure!(
                    had == &value,
                    "{:?} shard {:?} has rows with {key}={had} and {key}={value}; \
                     one row per run per shard cannot hold both",
                    self.what,
                    self.shard,
                );
            }
            self.params.insert(key, value);
        }

        let wall = row.wall_s().unwrap_or(0.0);
        let at = self.steps.entry(step).or_default();
        match row.batched() {
            true => at.1 = at.1.max(wall),
            false => at.0 += wall,
        }

        Ok(())
    }

    fn row(&self) -> Row {
        let wall = self
            .steps
            .values()
            .map(|(serial, batch)| serial + batch)
            .sum();

        Row {
            name: match self.stage {
                true => String::new(),
                false => self.what.clone(),
            },
            tool: self.tool.clone(),
            shard: self.shard.clone(),
            stage: match self.stage {
                true => self.what.clone(),
                false => String::new(),
            },
            params: self.params.clone(),
            wall_s: Some(wall),
        }
    }
}

/// The file [`Ledger::write`] writes, in a pipeline directory.
pub fn path(dir: &Path) -> PathBuf {
    dir.join(FILE)
}

/// Drops the ledger of whatever ran here last, before a pipeline replaces the
/// results it describes.
///
/// A pipeline that then dies leaves no ledger, and [`Ledger::load`] falls back
/// to the manifest of the run that just failed rather than reading the one
/// before it.
pub fn clear(dir: &Path) {
    std::fs::remove_file(path(dir)).ok();
}

/// Writes the ledger of a pipeline that has just finished.
pub fn record(dir: &Path) -> anyhow::Result<()> {
    let out = path(dir);
    Ledger::distill(&dir.join("manifest.tbl"))?.write(&out)?;

    println!("wrote {}", out.display());
    Ok(())
}

/// What a caller should say about the runs a pipeline did not finish.
pub fn warn(failed: &[(String, String)], what: &str) {
    if failed.is_empty() {
        return;
    }

    let listed: Vec<String> = failed
        .iter()
        .map(|(name, shard)| match shard.is_empty() {
            true => name.clone(),
            false => format!("{name}.{shard}"),
        })
        .collect();

    eprintln!(
        "warning: leaving out {} {what} that did not finish: {}",
        listed.len(),
        listed.join(", ")
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    // one directory per test: they run in parallel, and a shared one gets
    // removed out from under its neighbours
    fn distilled(name: &str, text: &str) -> anyhow::Result<Ledger> {
        let dir = std::env::temp_dir().join(format!("util-ledger-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("manifest.tbl");
        std::fs::write(&path, text).unwrap();

        Ledger::distill(&path)
    }

    fn manifest(name: &str, text: &str) -> Ledger {
        distilled(name, text).unwrap()
    }

    #[test]
    fn a_batched_step_is_its_longest_command_and_steps_add() {
        // one hmmer run over one shard: two parts alongside each other, then
        // the cat that joins them
        let runs = manifest(
            "batched",
            "# step cmd   name  shard tool  wall(s) exit\n\
             # ---- ----- ----- ----- ----- ------- ----\n\
             [1]    -     -     -     -     0.55    -\n\
             \x20  ||    hmmer hmmer 1     hmmer 0.25    0\n\
             \x20  ||    hmmer hmmer 1     hmmer 0.30    0\n\
             [2]    cat   hmmer 1     hmmer 0.05    0\n",
        );

        let columns = runs.columns().unwrap();
        assert_eq!(columns.len(), 1);
        assert_eq!(columns[0].wall_s, 0.35);
    }

    #[test]
    fn two_batched_steps_in_one_shard_both_count() {
        // what `Wall` could not express: a bucket holding more than one
        // batched step reported the longer of the two
        let runs = manifest(
            "two-steps",
            "# step cmd  name shard tool wall(s) exit\n\
             # ---- ---- ---- ----- ---- ------- ----\n\
             [1]    -    -    -     -    0.30    -\n\
             \x20  ||   nail a    1     nail 0.30    0\n\
             [2]    -    -    -     -    0.40    -\n\
             \x20  ||   nail a    1     nail 0.40    0\n",
        );

        assert_eq!(runs.columns().unwrap()[0].wall_s, 0.70);
    }

    #[test]
    fn commands_sharing_a_name_with_no_shard_axis_are_one_run() {
        // pid's shape: mmseqs' search and its conversion are one run, and the
        // benchmark has no shard axis at all
        let runs = manifest(
            "shardless",
            "# step cmd         name       mode tool   wall(s) exit\n\
             # ---- ----------- ---------- ---- ------ ------- ----\n\
             [1]    mmseqs      mmseqs.prf prf  mmseqs 4.00    0\n\
             [2]    convertalis mmseqs.prf prf  mmseqs 1.00    0\n",
        );

        let columns = runs.columns().unwrap();
        assert_eq!(columns.len(), 1);
        assert_eq!(columns[0].wall_s, 5.00);
        assert_eq!(columns[0].params.get("mode").unwrap(), "prf");
        assert_eq!(runs.shards(), vec![String::new()]);
    }

    #[test]
    fn a_stage_is_kept_apart_from_the_runs() {
        let runs = manifest(
            "stage",
            "# step cmd  name shard stage tool wall(s) exit\n\
             # ---- ---- ---- ----- ----- ---- ------- ----\n\
             [1]    seed -    1     seed  -    2.16    0\n\
             [2]    nail a    1     -     nail 0.12    0\n",
        );

        assert_eq!(runs.runs().count(), 1);
        let seeds: Vec<f64> = runs.stage("seed").filter_map(|row| row.wall_s).collect();
        assert_eq!(seeds, vec![2.16]);
    }

    #[test]
    fn a_failed_command_is_left_out_and_named() {
        let runs = manifest(
            "failed",
            "# step cmd  name shard tool wall(s) exit\n\
             # ---- ---- ---- ----- ---- ------- ----\n\
             [1]    nail a    1     nail 0.12    0\n\
             [2]    nail a    2     nail 0.00    1\n",
        );

        assert_eq!(runs.shards(), vec!["1".to_string()]);
        assert_eq!(runs.failed(), [("a".to_string(), "2".to_string())]);
    }

    #[test]
    fn rows_of_one_shard_that_disagree_on_a_param_are_refused() {
        let err = distilled(
            "reps",
            "# step cmd  name shard rep tool wall(s) exit\n\
             # ---- ---- ---- ----- --- ---- ------- ----\n\
             [1]    nail a    1     1   nail 0.12    0\n\
             [2]    nail a    1     2   nail 0.11    0\n",
        )
        .err()
        .unwrap()
        .to_string();

        assert!(err.contains("rep=1"), "{err}");
        assert!(err.contains("rep=2"), "{err}");
    }

    #[test]
    fn a_whole_run_figure_replaces_the_sum_over_shards() {
        let rows = vec![
            Row {
                name: "hmmer".to_string(),
                tool: "hmmer".to_string(),
                shard: "1".to_string(),
                stage: String::new(),
                params: BTreeMap::new(),
                wall_s: None,
            },
            Row {
                name: "hmmer".to_string(),
                tool: "hmmer".to_string(),
                shard: "2".to_string(),
                stage: String::new(),
                params: BTreeMap::new(),
                wall_s: None,
            },
            Row {
                name: "hmmer".to_string(),
                tool: "hmmer".to_string(),
                shard: EVERY_SHARD.to_string(),
                stage: String::new(),
                params: BTreeMap::new(),
                wall_s: Some(2043.77),
            },
        ];

        let columns = Ledger::from_rows(rows).columns().unwrap();
        assert_eq!(columns[0].wall_s, 2043.77);
        assert_eq!(columns[0].shards, vec!["1".to_string(), "2".to_string()]);
    }
}
