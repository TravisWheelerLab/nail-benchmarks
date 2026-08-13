//! A sink that writes the summary table.
//!
//! Two kinds of row, told apart by the first column. A step row carries the
//! step's name and its wall clock; the command rows under it carry a `|` or `||`
//! instead of a name, and their own numbers. A step of one command collapses to
//! a single row, since the two would otherwise say the same thing twice.
//!
//! Rows are buffered per step, because a column's width is not known until the
//! last cell in it has arrived. What happens at the end of a step depends on the
//! mode: [`Table::new`] writes that step's block right away, [`Table::whole`]
//! holds everything back so the entire table lines up as one.

use std::collections::BTreeSet;
use std::path::PathBuf;

use anyhow::Context;

use crate::cmd::Cmd;
use crate::execute::{Status, Timing};
use crate::sink::Sink;
use crate::step::Step;

/// Marks a command that ran after the one above it.
const SERIAL: &str = "|";
/// Marks a command that ran alongside the others in its step.
const BATCH: &str = "||";

/// Writes a table of everything the pipeline ran.
#[derive(Debug)]
pub struct Table {
    path: PathBuf,
    whole: bool,
    /// Every label key and tag any command carries, worked out up front so all
    /// blocks end up with the same columns whatever order the labels turn up in.
    header: Vec<String>,
    dimensions: usize,
    keys: Vec<String>,
    tags: Vec<String>,
    /// The current step's command rows, waiting for its step row.
    pending: Vec<Vec<Cell>>,
    /// Blocks already written, or every row so far in whole mode.
    done: Vec<Vec<Cell>>,
    text: String,
}

impl Table {
    /// Write each step's block as that step finishes. Columns are the same
    /// throughout, but each block is padded to its own contents, so blocks do
    /// not line up with each other.
    pub fn new(path: impl Into<PathBuf>) -> Table {
        Table {
            path: path.into(),
            whole: false,
            header: Vec::new(),
            dimensions: 0,
            keys: Vec::new(),
            tags: Vec::new(),
            pending: Vec::new(),
            done: Vec::new(),
            text: String::new(),
        }
    }

    /// Hold everything until the run is over and write one table, every row
    /// padded to the same widths. Nothing reaches the file until the end.
    pub fn whole(path: impl Into<PathBuf>) -> Table {
        Table {
            whole: true,
            ..Table::new(path)
        }
    }

    fn flush(&mut self) -> anyhow::Result<()> {
        if let Some(dir) = self.path.parent() {
            std::fs::create_dir_all(dir)?;
        }
        std::fs::write(&self.path, &self.text)
            .with_context(|| format!("failed to write {}", self.path.display()))
    }
}

impl Sink for Table {
    fn start(&mut self, steps: &[Step]) -> anyhow::Result<()> {
        let cmds = || steps.iter().flat_map(|s| s.cmds());

        self.keys = cmds()
            .flat_map(|c| c.labels.keys().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        self.tags = cmds()
            .flat_map(|c| c.tags.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        self.dimensions = self.keys.len() + self.tags.len();

        self.header = vec!["step".into(), "cmd".into()];
        self.header.extend(self.keys.iter().cloned());
        self.header.extend(self.tags.iter().cloned());
        self.header.extend(
            [
                "wall(s)", "user(s)", "sys(s)", "max_rss", "exit", "status", "argv",
            ]
            .iter()
            .map(|s| s.to_string()),
        );

        self.text.clear();
        self.flush()
    }

    fn record(&mut self, step: &Step, cmd: &Cmd) -> anyhow::Result<()> {
        // a step of one would just repeat itself, so the step keeps the first
        // column and the command fills in the rest of the row
        let first = if step.cmds().len() == 1 {
            Cell::left(step.name())
        } else {
            Cell::right(if step.batched() { BATCH } else { SERIAL })
        };
        self.pending
            .push(cmd_row(first, cmd, &self.keys, &self.tags));
        Ok(())
    }

    fn step_done(&mut self, step: &Step) -> anyhow::Result<()> {
        let mut block: Vec<Vec<Cell>> = Vec::new();
        if step.cmds().len() != 1 {
            block.push(step_row(step, self.dimensions));
        }
        block.append(&mut self.pending);

        if self.whole {
            self.done.append(&mut block);
            return Ok(());
        }

        if !self.text.is_empty() {
            self.text.push('\n');
        }
        self.text.push_str(&render(&self.header, &block));
        self.flush()
    }

    fn finish(&mut self) -> anyhow::Result<()> {
        if self.whole {
            let rows = std::mem::take(&mut self.done);
            self.text = render(&self.header, &rows);
            self.flush()?;
        }
        Ok(())
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

    // exit and status belong to commands; a count or a rollup here would be a
    // different quantity sharing a column
    cells.extend(std::iter::repeat_n(Cell::left("-"), 3));
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
            cells.push(Cell::left("ok"));
        }
        // nothing was measured either way; status is what separates "could not
        // start" from "never got its turn", which leaves exit to hold exit codes
        // and nothing else
        Status::Failed(_) => {
            cells.extend(std::iter::repeat_n(Cell::left("-"), 5));
            cells.push(Cell::left("err"));
        }
        Status::Skipped => {
            cells.extend(std::iter::repeat_n(Cell::left("-"), 5));
            cells.push(Cell::left("skip"));
        }
        Status::NotRun => cells.extend(std::iter::repeat_n(Cell::left("-"), 6)),
    }

    cells.push(Cell::left(cmd.line()));
    cells
}

/// A cell and which side its padding goes on.
#[derive(Clone, Debug)]
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
