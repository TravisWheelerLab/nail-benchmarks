//! The runs table: one row per (run, search), written as a benchmark executes
//! and read back by its analysis.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use anyhow::{bail, Context};

use crate::config::Run;
use crate::exec::Timing;

pub const FILE_NAME: &str = "runs.tbl";

/// How often the table is rewritten while a run is in progress.
///
/// Column widths are not known until the last row, so every rewrite emits the
/// whole file. Doing that per row is quadratic in bytes written. Rewriting on a
/// timer instead bounds the loss from a kill to the last few seconds.
const FLUSH_EVERY: Duration = Duration::from_secs(5);

// ------------------------------------------------------------------ writing

/// Accumulates rows and rewrites the table as it goes.
///
/// Columns are space-padded to the widest cell and the header is commented out
/// with a dashed separator, matching the layout nail uses for its own `.tbl`
/// output. Widths are not known until the last row arrives, so rows are kept in
/// memory and the whole file is rewritten, leaving the file on disk complete
/// and aligned even if a long run is interrupted.
pub struct Writer {
    path: PathBuf,
    sweep_columns: Vec<String>,
    header: Vec<String>,
    rows: Vec<Vec<String>>,
    last_flush: Instant,
}

impl Writer {
    /// Create the table at an explicit path.
    ///
    /// It does not have to live beside the hit tables: putting it above them
    /// keeps it from being wiped when a run clears results.
    pub fn create(path: impl AsRef<Path>, sweep_columns: Vec<String>) -> anyhow::Result<Self> {
        let path = path.as_ref().to_path_buf();
        if let Some(dir) = path.parent() {
            std::fs::create_dir_all(dir)?;
        }

        let mut header: Vec<String> = ["name", "tool", "query"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        header.extend(sweep_columns.iter().cloned());
        header.extend(
            ["threads", "search", "wall(s)", "user(s)", "sys(s)", "max_rss", "exit", "cmd"]
                .iter()
                .map(|s| s.to_string()),
        );

        let mut table = Writer {
            path,
            sweep_columns,
            header,
            rows: Vec::new(),
            last_flush: Instant::now(),
        };

        table.flush()?;
        Ok(table)
    }

    pub fn append(
        &mut self,
        run: &Run,
        search: &str,
        timing: &Timing,
        cmd: &str,
    ) -> anyhow::Result<()> {
        let mut row: Vec<String> = vec![
            run.name.clone(),
            run.tool.clone(),
            run.var_str("query").unwrap_or_else(|| "-".to_string()),
        ];

        // a run only carries the axes its own block swept, so anything from
        // another block's union shows up as "-"
        for col in &self.sweep_columns {
            row.push(
                run.var(col)
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
            );
        }

        row.push(run.threads.to_string());
        row.push(if search.is_empty() { "-" } else { search }.to_string());
        row.push(format!("{:.2}", timing.wall_s));
        row.push(format!("{:.2}", timing.user_s));
        row.push(format!("{:.2}", timing.sys_s));
        row.push(format_bytes(timing.max_rss_kb));
        row.push(timing.exit.to_string());
        row.push(cmd.to_string());

        self.rows.push(row);

        if self.last_flush.elapsed() >= FLUSH_EVERY {
            self.flush()?;
        }

        Ok(())
    }

    /// Write the table out unconditionally. Call this when the run ends, so the
    /// file reflects every row rather than the last timed rewrite.
    pub fn flush(&mut self) -> anyhow::Result<()> {
        std::fs::write(&self.path, render(&self.header, &self.rows))
            .with_context(|| format!("failed to write {}", self.path.display()))?;

        self.last_flush = Instant::now();
        Ok(())
    }

    pub fn path(&self) -> &Path {
        &self.path
    }
}

/// Lay out a commented header, a dashed separator, and the data rows, with
/// every column padded to its widest cell.
fn render(header: &[String], rows: &[Vec<String>]) -> String {
    // the "# " comment marker is absorbed into the first column's width so the
    // header labels stay lined up over the data below them
    let mut head = header.to_vec();
    head[0] = format!("# {}", head[0]);

    let mut widths: Vec<usize> = head.iter().map(|h| h.chars().count()).collect();
    for row in rows {
        for (i, cell) in row.iter().enumerate() {
            widths[i] = widths[i].max(cell.chars().count());
        }
    }

    let last = widths.len() - 1;
    let mut sep: Vec<String> = widths.iter().map(|w| "-".repeat(*w)).collect();
    sep[0] = format!("# {}", "-".repeat(widths[0].saturating_sub(2)));
    // cmd is unpadded and hundreds of characters wide; underline the label
    // rather than the whole column
    sep[last] = "-".repeat(header[last].chars().count());

    let mut out = String::new();
    write_row(&mut out, &head, &widths);
    write_row(&mut out, &sep, &widths);
    for row in rows {
        write_row(&mut out, row, &widths);
    }

    out
}

/// Peak memory in binary units, kept to about three significant figures so the
/// column reads at a glance: 940KiB, 10.4MiB, 1.02GiB.
fn format_bytes(kib: i64) -> String {
    if kib < 0 {
        return "-".to_string();
    }

    const STEP: f64 = 1024.0;
    let (value, unit) = match kib as f64 {
        v if v < STEP => (v, "KiB"),
        v if v < STEP * STEP => (v / STEP, "MiB"),
        v => (v / (STEP * STEP), "GiB"),
    };

    if value >= 100.0 {
        format!("{value:.0}{unit}")
    } else if value >= 10.0 {
        format!("{value:.1}{unit}")
    } else {
        format!("{value:.2}{unit}")
    }
}

fn write_row(out: &mut String, cells: &[String], widths: &[usize]) {
    let last = cells.len().saturating_sub(1);
    for (i, cell) in cells.iter().enumerate() {
        if i == last {
            // the final column is cmd, which is long and variable; padding it
            // would only add trailing whitespace
            out.push_str(cell);
        } else {
            let pad = widths[i].saturating_sub(cell.chars().count());
            out.push_str(cell);
            out.extend(std::iter::repeat_n(' ', pad + 1));
        }
    }
    out.push('\n');
}

// ------------------------------------------------------------------ reading

/// One executed (run, search) pair.
#[derive(Clone, Debug)]
pub struct Row {
    pub name: String,
    pub tool: String,
    /// The search's label, or empty when the benchmark had only one.
    pub search: String,
    pub wall_s: f32,
    pub exit: i32,
}

/// A runs table read back, with the hit tables it indexes.
#[derive(Clone, Debug)]
pub struct Runs {
    dir: PathBuf,
    rows: Vec<Row>,
}

impl Runs {
    /// Read `<dir>/runs.tbl`, with the hit tables alongside it.
    pub fn from_dir(dir: impl AsRef<Path>) -> anyhow::Result<Self> {
        let dir = dir.as_ref();
        Self::load(dir.join(FILE_NAME), dir)
    }

    /// Read a runs table that does not sit beside the hit tables it indexes.
    pub fn load(table: impl AsRef<Path>, results_dir: impl AsRef<Path>) -> anyhow::Result<Self> {
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

            let idx = |name: &str| -> anyhow::Result<usize> {
                columns
                    .iter()
                    .position(|c| c == name)
                    .with_context(|| format!("{} has no `{name}` column", path.display()))
            };

            // cmd is last and contains spaces, but every column read here comes
            // before it, so positional splitting is safe
            let f: Vec<&str> = line.split_whitespace().collect();
            let get = |i: usize| -> anyhow::Result<&str> {
                f.get(i).copied().context("short row in runs table")
            };

            let search = get(idx("search")?)?;
            rows.push(Row {
                name: get(idx("name")?)?.to_string(),
                tool: get(idx("tool")?)?.to_string(),
                search: if search == "-" { String::new() } else { search.to_string() },
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
    pub fn only_for_tool(&self, tool: &str) -> anyhow::Result<&str> {
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

    /// Search labels a given run was executed against, in order.
    pub fn searches(&self, name: &str) -> Vec<&str> {
        self.rows
            .iter()
            .filter(|r| r.name == name)
            .map(|r| r.search.as_str())
            .collect()
    }

    /// Labels common to every named run, so callers comparing tools only walk
    /// pairs that actually exist.
    pub fn shared_searches(&self, names: &[&str]) -> Vec<String> {
        let mut common: Option<BTreeSet<String>> = None;

        for name in names {
            let set: BTreeSet<String> = self
                .rows
                .iter()
                .filter(|r| &r.name == name && r.exit == 0)
                .map(|r| r.search.clone())
                .collect();

            common = Some(match common {
                None => set,
                Some(prev) => prev.intersection(&set).cloned().collect(),
            });
        }

        let mut out: Vec<String> = common.unwrap_or_default().into_iter().collect();
        // labels are usually numeric shard ids, which sort badly as strings
        out.sort_by(|a, b| match (a.parse::<u64>(), b.parse::<u64>()) {
            (Ok(x), Ok(y)) => x.cmp(&y),
            _ => a.cmp(b),
        });
        out
    }

    /// Where a run's hit table for one search lives.
    pub fn table_path(&self, name: &str, search: &str) -> PathBuf {
        self.output_path(name, search, "tbl")
    }

    /// Outputs are named for the run and the search that produced them, which
    /// is what [`crate::Paths::out`] writes.
    pub fn output_path(&self, name: &str, search: &str, ext: &str) -> PathBuf {
        if search.is_empty() {
            self.dir.join(format!("{name}.{ext}"))
        } else {
            self.dir.join(format!("{name}.{search}.{ext}"))
        }
    }

    pub fn wall_s(&self, name: &str, search: &str) -> Option<f32> {
        self.rows
            .iter()
            .find(|r| r.name == name && r.search == search)
            .map(|r| r.wall_s)
    }

    /// Mean wall-clock seconds across every search a run covered.
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

    fn strings(cells: &[&str]) -> Vec<String> {
        cells.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn memory_reads_in_binary_units() {
        assert_eq!(format_bytes(940), "940KiB");
        assert_eq!(format_bytes(10_649), "10.4MiB");
        assert_eq!(format_bytes(107_800), "105MiB");
        assert_eq!(format_bytes(1_066_344), "1.02GiB");
        // exactly on a boundary should step up rather than read as 1024KiB
        assert_eq!(format_bytes(1024), "1.00MiB");
    }

    #[test]
    fn rows_pad_to_the_widest_cell_in_each_column() {
        let header = strings(&["name", "s", "cmd"]);
        let rows = vec![
            strings(&["a", "12.0", "/bin/x --y"]),
            strings(&["longer-name", "5.7", "/bin/z"]),
        ];

        let out = render(&header, &rows);
        let lines: Vec<&str> = out.lines().collect();

        assert_eq!(lines[0], "# name      s    cmd");
        assert_eq!(lines[1], "# --------- ---- ---");
        assert_eq!(lines[2], "a           12.0 /bin/x --y");
        assert_eq!(lines[3], "longer-name 5.7  /bin/z");

        // no line carries trailing padding
        for line in &lines {
            assert_eq!(line.trim_end(), *line, "trailing whitespace in {line:?}");
        }
    }

    #[test]
    fn header_marker_does_not_shift_the_columns() {
        let header = strings(&["name", "tool", "cmd"]);
        let rows = vec![
            strings(&["nail-s12.0.prf", "nail", "/bin/nail"]),
            strings(&["hmmer.seq", "hmmer", "/bin/hmmsearch"]),
        ];

        let out = render(&header, &rows);

        // the "# " marker eats into the first column rather than pushing
        // everything right, so the tool column starts at one offset on every
        // line: header, separator, and data alike
        let tool_starts: Vec<usize> = out
            .lines()
            .map(|line| {
                // skip the marker so header and data lines are measured alike
                let from = if line.starts_with("# ") { 2 } else { 0 };
                let gap = line[from..].find(' ').unwrap() + from;
                line[gap..].find(|c: char| c != ' ').unwrap() + gap
            })
            .collect();

        assert!(
            tool_starts.windows(2).all(|w| w[0] == w[1]),
            "columns not aligned across lines: {tool_starts:?}\n{out}"
        );
    }

    // one directory per test: they run in parallel, and a shared one gets
    // deleted out from under its neighbours
    fn write(name: &str, body: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("bm-table-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join(FILE_NAME), body).unwrap();
        dir
    }

    const BODY: &str = "\
# name                  tool   query s    threads search wall(s) user(s) sys(s) max_rss exit cmd
# --------------------- ------ ----- ---- ------- ------ ------ ------ ----- ---------- ---- ---
nail-s12.0-prog.prf     nail   prf   12.0 24      1      96.73  2068.7 0.5   1.02GiB    0    /bin/nail search --a b
nail-s12.0-prog.prf     nail   prf   12.0 24      2      95.83  2058.6 0.5   1.03GiB    0    /bin/nail search --a b
mmseqs-s7.5-ms2000.prf  mmseqs prf   7.5  24      1      7.01   161.0  0.4   636MiB     0    /bin/mmseqs search x y
mmseqs-s7.5-ms2000.prf  mmseqs prf   7.5  24      2      7.03   161.6  0.4   635MiB     0    /bin/mmseqs search x y
";

    #[test]
    fn names_survive_dots_in_parameters() {
        let dir = write("names", BODY);
        let runs = Runs::from_dir(&dir).unwrap();

        // splitting these on '.' would mangle s12.0
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
            runs.table_path("nail-s12.0-prog.prf", "2"),
            dir.join("nail-s12.0-prog.prf.2.tbl")
        );
        assert_eq!(runs.wall_s("nail-s12.0-prog.prf", "1"), Some(96.73));
        assert_eq!(
            runs.mean_wall_s("mmseqs-s7.5-ms2000.prf").unwrap(),
            (7.01 + 7.03) / 2.0
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn an_unlabelled_search_leaves_its_token_out_of_the_path() {
        let dir = write(
            "single",
            "\
# name      tool query threads search wall(s) user(s) sys(s) max_rss exit cmd
# --------- ---- ----- ------- ------ ------- ------- ------ ------- ---- ---
nail.prf    nail prf   24      -      1.00    2.00    0.10   1MiB    0    /bin/nail x
",
        );
        let runs = Runs::from_dir(&dir).unwrap();

        assert_eq!(runs.table_path("nail.prf", ""), dir.join("nail.prf.tbl"));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn shared_searches_are_the_intersection_in_numeric_order() {
        let dir = write("shared", BODY);
        let runs = Runs::from_dir(&dir).unwrap();

        let shared = runs.shared_searches(&["nail-s12.0-prog.prf", "mmseqs-s7.5-ms2000.prf"]);
        assert_eq!(shared, vec!["1", "2"]);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn ambiguous_tool_selection_is_rejected() {
        let dir = write("ambiguous", &format!(
            "{BODY}nail-s9.0-prog.prf       nail   prf   9.0  24      1      12.0   30.0   0.2   100KiB     0    /bin/nail x y\n"
        ));
        let runs = Runs::from_dir(&dir).unwrap();

        let err = runs.only_for_tool("nail").unwrap_err().to_string();
        assert!(err.contains("2 runs for tool"), "unexpected: {err}");
        assert_eq!(runs.only_for_tool("mmseqs").unwrap(), "mmseqs-s7.5-ms2000.prf");
        std::fs::remove_dir_all(&dir).ok();
    }
}
