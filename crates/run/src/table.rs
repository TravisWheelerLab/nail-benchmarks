//! Reading back the runs table this crate writes.
//!
//! Analysis code used to rediscover what had been run by globbing result files
//! and picking filenames apart. That is fragile — run names embed floats, so a
//! parameter like `s12.0` puts a dot in the middle of a name that was being
//! split on dots. The runs table already records every name, target, and
//! timing, so it is the authoritative index of a results directory.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

pub const FILE_NAME: &str = "runs.tbl";

/// One executed (run, target) pair.
#[derive(Clone, Debug)]
pub struct Row {
    pub name: String,
    pub tool: String,
    pub target: String,
    pub wall_s: f32,
    pub exit: i32,
}

#[derive(Clone, Debug)]
pub struct Runs {
    dir: PathBuf,
    rows: Vec<Row>,
}

impl Runs {
    /// Read `<dir>/runs.tbl`, with the hit tables alongside it.
    pub fn from_dir(dir: impl AsRef<Path>) -> Result<Self> {
        let dir = dir.as_ref();
        Self::load(dir.join(FILE_NAME), dir)
    }

    /// Read a runs table that does not sit beside the hit tables it indexes.
    pub fn load(table: impl AsRef<Path>, results_dir: impl AsRef<Path>) -> Result<Self> {
        let dir = results_dir.as_ref().to_path_buf();
        let path = table.as_ref().to_path_buf();
        let text = std::fs::read_to_string(&path).with_context(|| {
            format!(
                "failed to read {}; has this benchmark been run?",
                path.display()
            )
        })?;

        let mut columns: Vec<String> = Vec::new();
        let mut rows = Vec::new();

        for line in text.lines() {
            if let Some(rest) = line.strip_prefix("# ") {
                // the header names the columns; the separator row after it is
                // only dashes
                if columns.is_empty() && !rest.starts_with('-') {
                    columns = rest.split_whitespace().map(str::to_string).collect();
                }
                continue;
            }

            if line.trim().is_empty() {
                continue;
            }

            let idx = |name: &str| -> Result<usize> {
                columns
                    .iter()
                    .position(|c| c == name)
                    .with_context(|| format!("{} has no `{name}` column", path.display()))
            };

            // cmd is last and contains spaces, but every column read here comes
            // before it, so positional splitting is safe
            let f: Vec<&str> = line.split_whitespace().collect();
            let get = |i: usize| -> Result<&str> {
                f.get(i).copied().context("short row in runs table")
            };

            rows.push(Row {
                name: get(idx("name")?)?.to_string(),
                tool: get(idx("tool")?)?.to_string(),
                target: get(idx("target")?)?.to_string(),
                wall_s: get(idx("wall(s)")?)?.parse().unwrap_or(f32::NAN),
                exit: get(idx("exit")?)?.parse().unwrap_or(-1),
            });
        }

        if rows.is_empty() {
            bail!("no runs recorded in {}", path.display());
        }

        Ok(Runs { dir, rows })
    }

    pub fn rows(&self) -> &[Row] {
        &self.rows
    }

    /// Distinct run names, in the order they were executed.
    pub fn names(&self) -> Vec<&str> {
        let mut seen = BTreeSet::new();
        let mut out = Vec::new();
        for row in &self.rows {
            if seen.insert(row.name.as_str()) {
                out.push(row.name.as_str());
            }
        }
        out
    }

    /// The single run name belonging to `tool`, erroring if the config swept
    /// more than one so the caller has to say which it means.
    pub fn only_for_tool(&self, tool: &str) -> Result<&str> {
        let names: Vec<&str> = {
            let mut seen = BTreeSet::new();
            self.rows
                .iter()
                .filter(|r| r.tool == tool)
                .map(|r| r.name.as_str())
                .filter(|n| seen.insert(*n))
                .collect()
        };

        match names.len() {
            0 => bail!("no runs for tool {tool:?} in {}", self.dir.display()),
            1 => Ok(names[0]),
            _ => bail!(
                "{} runs for tool {tool:?}; pass one of: {}",
                names.len(),
                names.join(", ")
            ),
        }
    }

    /// Targets a given run was executed against, in order.
    pub fn targets(&self, name: &str) -> Vec<&str> {
        self.rows
            .iter()
            .filter(|r| r.name == name)
            .map(|r| r.target.as_str())
            .collect()
    }

    /// Targets common to every named run, so callers comparing tools only walk
    /// pairs that actually exist.
    pub fn shared_targets(&self, names: &[&str]) -> Vec<String> {
        let mut common: Option<BTreeSet<String>> = None;

        for name in names {
            let set: BTreeSet<String> = self
                .rows
                .iter()
                .filter(|r| &r.name == name && r.exit == 0)
                .map(|r| r.target.clone())
                .collect();

            common = Some(match common {
                None => set,
                Some(prev) => prev.intersection(&set).cloned().collect(),
            });
        }

        let mut out: Vec<String> = common.unwrap_or_default().into_iter().collect();
        // targets are usually numeric shard ids, which sort badly as strings
        out.sort_by(|a, b| match (a.parse::<u64>(), b.parse::<u64>()) {
            (Ok(x), Ok(y)) => x.cmp(&y),
            _ => a.cmp(b),
        });
        out
    }

    /// Where a run's hit table for one target lives, matching how the runner
    /// names its outputs.
    ///
    /// The `target` column holds a file name; outputs are keyed by its stem, so
    /// a row and the file it refers to always agree.
    pub fn table_path(&self, name: &str, target: &str) -> PathBuf {
        self.output_path(name, target, "tbl")
    }

    pub fn output_path(&self, name: &str, target: &str, ext: &str) -> PathBuf {
        if self.targets(name).len() <= 1 {
            return self.dir.join(format!("{name}.{ext}"));
        }

        let key = Path::new(target)
            .file_stem()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| target.to_string());

        self.dir.join(format!("{name}.{key}.{ext}"))
    }

    pub fn wall_s(&self, name: &str, target: &str) -> Option<f32> {
        self.rows
            .iter()
            .find(|r| r.name == name && r.target == target)
            .map(|r| r.wall_s)
    }

    /// Mean wall-clock seconds across every target a run covered.
    pub fn mean_wall_s(&self, name: &str) -> Option<f32> {
        let times: Vec<f32> = self
            .rows
            .iter()
            .filter(|r| r.name == name && r.exit == 0)
            .map(|r| r.wall_s)
            .collect();

        if times.is_empty() {
            return None;
        }
        Some(times.iter().sum::<f32>() / times.len() as f32)
    }

    pub fn dir(&self) -> &Path {
        &self.dir
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // one directory per test: they run in parallel, and a shared one gets
    // deleted out from under its neighbours
    fn write(name: &str, body: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("bm-table-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join(FILE_NAME), body).unwrap();
        dir
    }

    const BODY: &str = "\
# name                  tool   query s    threads target wall(s) user(s) sys(s) max_rss exit cmd
# --------------------- ------ ----- ---- ------- ------ ------ ------ ----- ---------- ---- ---
nail-s12.0-prog.prf     nail   prf   12.0 24      1.fa   96.73  2068.7 0.5   1.02GiB    0    /bin/nail search --a b
nail-s12.0-prog.prf     nail   prf   12.0 24      2.fa   95.83  2058.6 0.5   1.03GiB    0    /bin/nail search --a b
mmseqs-s7.5-ms2000.prf  mmseqs prf   7.5  24      1.fa   7.01   161.0  0.4   636MiB     0    /bin/mmseqs search x y
mmseqs-s7.5-ms2000.prf  mmseqs prf   7.5  24      2.fa   7.03   161.6  0.4   635MiB     0    /bin/mmseqs search x y
";

    #[test]
    fn names_survive_dots_in_parameters() {
        let dir = write("names", BODY);
        let runs = Runs::from_dir(&dir).unwrap();

        // the old code split filenames on '.', which mangled s12.0
        assert_eq!(
            runs.names(),
            vec!["nail-s12.0-prog.prf", "mmseqs-s7.5-ms2000.prf"]
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn table_paths_and_timings_resolve() {
        let dir = write("paths", BODY);
        let runs = Runs::from_dir(&dir).unwrap();

        assert_eq!(
            runs.table_path("nail-s12.0-prog.prf", "2.fa"),
            dir.join("nail-s12.0-prog.prf.2.tbl")
        );
        assert_eq!(runs.wall_s("nail-s12.0-prog.prf", "1.fa"), Some(96.73));
        assert_eq!(
            runs.mean_wall_s("mmseqs-s7.5-ms2000.prf").unwrap(),
            (7.01 + 7.03) / 2.0
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn shared_targets_are_the_intersection_in_numeric_order() {
        let dir = write("shared", BODY);
        let runs = Runs::from_dir(&dir).unwrap();

        let shared = runs.shared_targets(&["nail-s12.0-prog.prf", "mmseqs-s7.5-ms2000.prf"]);
        assert_eq!(shared, vec!["1.fa", "2.fa"]);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn ambiguous_tool_selection_is_rejected() {
        let dir = write("ambiguous", &format!(
            "{BODY}nail-s9.0-prog.prf       nail   prf   9.0  24      1.fa   12.0   30.0   0.2   100KiB     0    /bin/nail x y\n"
        ));
        let runs = Runs::from_dir(&dir).unwrap();

        let err = runs.only_for_tool("nail").unwrap_err().to_string();
        assert!(err.contains("2 runs for tool"), "unexpected: {err}");
        assert_eq!(runs.only_for_tool("mmseqs").unwrap(), "mmseqs-s7.5-ms2000.prf");
        std::fs::remove_dir_all(&dir).ok();
    }
}
