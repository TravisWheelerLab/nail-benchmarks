//! The summary table, written as the pipeline runs.
//!
//! Two kinds of row, told apart by the first column. A step row carries the
//! step's name and its wall clock; the command rows under it carry a `|` or `||`
//! instead of a name, and their own numbers. A step of one command collapses to
//! a single row, since the two would otherwise say the same thing twice.

use std::collections::BTreeSet;
use std::path::Path;

use anyhow::Context;

use crate::Step;
use crate::cmd::Cmd;
use crate::execute::{Status, Timing};

/// Marks a command that ran after the one above it.
const SERIAL: &str = "|";
/// Marks a command that ran alongside the others in its step.
const BATCH: &str = "||";

/// What a pipeline did: the steps it was given, with a [`crate::Status`] on
/// every command.
#[derive(Clone, Debug, Default)]
pub struct Report {
    pub(crate) steps: Vec<Step>,
}

impl Report {
    /// Every command, in the order they were meant to run, paired with the name
    /// of the step it belongs to.
    pub fn cmds(&self) -> impl Iterator<Item = (&str, &Cmd)> {
        self.steps
            .iter()
            .flat_map(|s| s.cmds().iter().map(move |c| (s.name(), c)))
    }

    /// How many commands could not be run or came back nonzero. Commands that
    /// never got their turn do not count.
    pub fn failed(&self) -> usize {
        self.cmds().filter(|(_, c)| c.status.failed()).count()
    }

    /// How long the pipeline spent inside steps, on the wall clock.
    pub fn wall_s(&self) -> f64 {
        self.steps.iter().filter_map(|s| s.wall_s()).sum()
    }

    /// The table, laid out with a commented header and padded columns.
    ///
    /// Label columns are the union of every key any command carries, so
    /// commands that did not set a key show a dash. Tag columns follow, reading
    /// `x` or `-`. `argv` is last and left unpadded, since it can run to
    /// hundreds of characters.
    pub fn render(&self) -> String {
        let keys: Vec<String> = self
            .cmds()
            .flat_map(|(_, c)| c.labels.keys().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();

        let tags: Vec<String> = self
            .cmds()
            .flat_map(|(_, c)| c.tags.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();

        let mut header: Vec<String> = vec!["step".into(), "cmd".into()];
        header.extend(keys.iter().cloned());
        header.extend(tags.iter().cloned());
        header.extend(
            ["wall(s)", "user(s)", "sys(s)", "max_rss", "exit", "argv"]
                .iter()
                .map(|s| s.to_string()),
        );

        let mut rows: Vec<Vec<Cell>> = Vec::new();
        for step in &self.steps {
            // a step of one would just repeat itself, so the step keeps the
            // first column and the command fills in the rest of the row
            if let [only] = step.cmds() {
                rows.push(cmd_row(Cell::left(step.name()), only, &keys, &tags));
                continue;
            }

            rows.push(step_row(step, keys.len() + tags.len()));
            let marker = if step.batched() { BATCH } else { SERIAL };
            for cmd in step.cmds() {
                rows.push(cmd_row(Cell::right(marker), cmd, &keys, &tags));
            }
        }

        render(&header, &rows)
    }

    pub fn write(&self, path: impl AsRef<Path>) -> anyhow::Result<()> {
        let path = path.as_ref();
        if let Some(dir) = path.parent() {
            std::fs::create_dir_all(dir)?;
        }
        std::fs::write(path, self.render())
            .with_context(|| format!("failed to write {}", path.display()))
    }

    pub fn print(&self) {
        print!("{}", self.render());
    }
}

/// The step's own line: measured wall clock, and its commands' CPU added up.
fn step_row(step: &Step, dimensions: usize) -> Vec<Cell> {
    let mut cells = vec![Cell::left(step.name()), Cell::left("-")];
    cells.extend(std::iter::repeat_n(Cell::left("-"), dimensions));

    cells.push(Cell::left(match step.wall_s() {
        Some(wall) => format!("{wall:.2}"),
        None => "-".to_string(),
    }));

    let timings: Vec<&Timing> = step
        .cmds()
        .iter()
        .filter_map(|c| match c.status() {
            Status::Finished(t) => Some(t),
            _ => None,
        })
        .collect();

    if timings.is_empty() {
        cells.extend(std::iter::repeat_n(Cell::left("-"), 3));
    } else {
        // summed CPU against measured wall is what shows whether a batch
        // actually bought anything
        cells.push(Cell::left(format!(
            "{:.2}",
            timings.iter().map(|t| t.user_s).sum::<f64>()
        )));
        cells.push(Cell::left(format!(
            "{:.2}",
            timings.iter().map(|t| t.sys_s).sum::<f64>()
        )));
        // the largest any one process got, which is not the same as the most
        // the step held at once — wait4 cannot tell us that
        cells.push(Cell::left(format_bytes(
            timings.iter().map(|t| t.max_rss_kb).max().unwrap_or(-1),
        )));
    }

    // exit codes belong to commands; a count here would be a different quantity
    // sharing a column
    cells.push(Cell::left("-"));
    cells.push(Cell::left("-"));
    cells
}

/// One command's line. `first` is the step name for a collapsed step of one,
/// and a right-aligned `|` or `||` otherwise.
fn cmd_row(first: Cell, cmd: &Cmd, keys: &[String], tags: &[String]) -> Vec<Cell> {
    let mut cells = vec![first, Cell::left(cmd.name())];
    cells.extend(
        keys.iter()
            .map(|k| Cell::left(cmd.labels.get(k).map(String::as_str).unwrap_or("-"))),
    );
    cells.extend(
        tags.iter()
            .map(|t| Cell::left(if cmd.tags.contains(t) { "x" } else { "-" })),
    );

    match &cmd.status {
        Status::Finished(t) => {
            cells.push(Cell::left(format!("{:.2}", t.wall_s)));
            cells.push(Cell::left(format!("{:.2}", t.user_s)));
            cells.push(Cell::left(format!("{:.2}", t.sys_s)));
            cells.push(Cell::left(format_bytes(t.max_rss_kb)));
            cells.push(Cell::left(t.exit.to_string()));
        }
        // nothing was measured either way; the exit cell is what separates
        // "could not start" from "never got its turn"
        Status::Failed(_) => {
            cells.extend(std::iter::repeat_n(Cell::left("-"), 4));
            cells.push(Cell::left("err"));
        }
        Status::NotRun => cells.extend(std::iter::repeat_n(Cell::left("-"), 5)),
    }

    cells.push(Cell::left(cmd.line()));
    cells
}

/// A cell and which side its padding goes on.
#[derive(Clone)]
struct Cell {
    text: String,
    right: bool,
}

impl Cell {
    fn left(text: impl Into<String>) -> Cell {
        Cell {
            text: text.into(),
            right: false,
        }
    }

    fn right(text: impl Into<String>) -> Cell {
        Cell {
            text: text.into(),
            right: true,
        }
    }
}

/// Lay out a commented header, a dashed separator, and the rows, every column
/// padded to its widest cell.
fn render(header: &[String], rows: &[Vec<Cell>]) -> String {
    // the "# " marker is absorbed into the first column's width, so the header
    // labels stay lined up over the data below them
    let mut head: Vec<Cell> = header.iter().map(Cell::left).collect();
    head[0] = Cell::left(format!("# {}", header[0]));

    let mut widths: Vec<usize> = head.iter().map(|c| c.text.chars().count()).collect();
    for row in rows {
        for (i, cell) in row.iter().enumerate() {
            widths[i] = widths[i].max(cell.text.chars().count());
        }
    }

    let last = widths.len() - 1;
    let mut sep: Vec<Cell> = widths.iter().map(|w| Cell::left("-".repeat(*w))).collect();
    sep[0] = Cell::left(format!("# {}", "-".repeat(widths[0].saturating_sub(2))));
    // argv is unpadded, so underline the label rather than the whole column
    sep[last] = Cell::left("-".repeat(header[last].chars().count()));

    let mut out = String::new();
    write_row(&mut out, &head, &widths);
    write_row(&mut out, &sep, &widths);
    for row in rows {
        write_row(&mut out, row, &widths);
    }
    out
}

fn write_row(out: &mut String, cells: &[Cell], widths: &[usize]) {
    let last = cells.len() - 1;
    for (i, cell) in cells.iter().enumerate() {
        if i == last {
            out.push_str(&cell.text);
            break;
        }
        let pad = " ".repeat(widths[i].saturating_sub(cell.text.chars().count()));
        if cell.right {
            out.push_str(&pad);
            out.push_str(&cell.text);
        } else {
            out.push_str(&cell.text);
            out.push_str(&pad);
        }
        out.push(' ');
    }
    out.push('\n');
}

/// Peak memory in binary units, to about three significant figures so the
/// column reads at a glance: 940KiB, 10.4MiB, 1.02GiB.
pub(crate) fn format_bytes(kib: i64) -> String {
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
