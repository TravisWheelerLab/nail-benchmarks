//! A sink that writes the summary table.
//!
//! Two kinds of row, told apart by the first column. A step row carries the
//! step's name and its wall clock; the command rows under it carry a `|` or `||`
//! instead of a name, and their own numbers. A step of one command collapses to
//! a single row, since the two would otherwise say the same thing twice.
//!
//! A block is built when its step finishes, because a column's width is not
//! known until the last cell in it has arrived. [`Mode`] decides what happens
//! to it then.

use std::collections::BTreeSet;
use std::path::PathBuf;

use anyhow::Context;

use crate::cmd::Cmd;
use crate::execute::{Status, Timing};
use crate::sink::Sink;
use crate::step::{Step, Strategy};

/// Marks a command that ran after the one above it.
const SERIAL: &str = "|";
/// Marks a command that ran alongside the others in its step.
const BATCH: &str = "||";

/// Whether every block gets its own header.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Headers {
    /// One header at the top of the file. Its widths are the floor for every
    /// block, so blocks line up with it and with each other until some value
    /// turns out wider than the label above it.
    Once,
    /// A header on every block, so each block reads on its own.
    Each,
}

/// How the table is laid out and when it reaches the file.
///
/// The combinations that make no sense cannot be written down: there is nothing
/// to decide about headers when the file holds one block, and ragged columns
/// force a header on every block.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    /// Hold everything back and write one table with every row padded alike.
    /// Nothing reaches the file until the run is over.
    Whole,
    /// Write each step's block as that step finishes. Every block carries the
    /// same columns, padded to its own contents. The default, with
    /// [`Headers::Once`].
    Blocks { headers: Headers },
    /// Write each step's block as that step finishes, each carrying only the
    /// columns its own commands use.
    Ragged,
}

/// Writes a table of everything the pipeline ran.
#[derive(Debug)]
pub struct Table {
    path: PathBuf,
    mode: Mode,
    /// Every field key and tag any command carries, worked out before anything
    /// runs so blocks agree on their columns whatever order fields turn up in.
    /// Unused by [`Mode::Ragged`], which asks each step instead.
    keys: Vec<String>,
    tags: Vec<String>,
    /// Every row so far, for [`Mode::Whole`].
    rows: Vec<Vec<Cell>>,
    /// Column widths every block starts from, worked out before the run from
    /// everything already known: names, fields, tags, argv. Only the numbers are
    /// missing, and their headings are wider than they usually are. Empty for
    /// the modes that do not share widths between blocks.
    floor: Vec<usize>,
    text: String,
}

impl Default for Mode {
    fn default() -> Mode {
        Mode::Blocks {
            headers: Headers::Once,
        }
    }
}

impl Table {
    pub fn new(path: impl Into<PathBuf>) -> Table {
        Table {
            path: path.into(),
            mode: Mode::default(),
            keys: Vec::new(),
            tags: Vec::new(),
            rows: Vec::new(),
            floor: Vec::new(),
            text: String::new(),
        }
    }

    pub fn mode(mut self, mode: Mode) -> Table {
        self.mode = mode;
        self
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
        let cmds: Vec<&Cmd> = steps.iter().flat_map(|s| s.cmds()).collect();
        (self.keys, self.tags) = dimensions(&cmds);

        let head = header(&self.keys, &self.tags);

        // Ragged blocks are meant to differ, and Whole renders in one go, so
        // neither has anything to share.
        if !matches!(self.mode, Mode::Ragged | Mode::Whole) {
            self.floor = vec![0; head.len()];
            let blocks: Vec<Vec<Cell>> = steps
                .iter()
                .flat_map(|s| block(s, &self.keys, &self.tags))
                .collect();
            widen(&mut self.floor, &head_cells(&head));
            for row in &blocks {
                widen(&mut self.floor, row);
            }
        }

        self.text.clear();
        if matches!(
            self.mode,
            Mode::Blocks {
                headers: Headers::Once
            }
        ) {
            self.text = render(&head, &[], true, &self.floor);
        }
        self.flush()
    }

    // no `record`: a step's rows are built from the step itself once it is done,
    // which keeps them in the order the commands were declared rather than the
    // order a batch happened to finish them in

    fn step_done(&mut self, step: &Step) -> anyhow::Result<()> {
        let (keys, tags) = match self.mode {
            Mode::Ragged => dimensions(&step.cmds().iter().collect::<Vec<_>>()),
            _ => (self.keys.clone(), self.tags.clone()),
        };

        let mut rows = block(step, &keys, &tags);

        let show_header = match self.mode {
            Mode::Whole => {
                self.rows.append(&mut rows);
                return Ok(());
            }
            Mode::Blocks {
                headers: Headers::Once,
            } => false,
            _ => true,
        };

        self.text.push_str(&render(
            &header(&keys, &tags),
            &rows,
            show_header,
            &self.floor,
        ));
        self.flush()
    }

    fn finish(&mut self) -> anyhow::Result<()> {
        if self.mode == Mode::Whole {
            let rows = std::mem::take(&mut self.rows);
            self.text = render(&header(&self.keys, &self.tags), &rows, true, &[]);
            self.flush()?;
        }
        Ok(())
    }
}

/// One step's rows: its own line, then a line per command.
fn block(step: &Step, keys: &[String], tags: &[String]) -> Vec<Vec<Cell>> {
    let mut rows = Vec::new();

    // a step of one would just repeat itself, so it gets no line of its own and
    // keeps the first column instead, with its command filling in the rest
    let first = if step.cmds().len() == 1 {
        Cell::left(step.label())
    } else {
        rows.push(step_row(step, keys.len() + tags.len()));
        Cell::right(match step.strategy() {
            Strategy::Serial => SERIAL,
            Strategy::Batch { .. } => BATCH,
        })
    };

    for cmd in step.cmds() {
        rows.push(cmd_row(first.clone(), cmd, keys, tags));
    }
    rows
}

/// The field keys and tags these commands carry, each in sorted order.
fn dimensions(cmds: &[&Cmd]) -> (Vec<String>, Vec<String>) {
    let keys = cmds
        .iter()
        .flat_map(|c| c.fields().keys().cloned())
        .collect::<BTreeSet<_>>();
    let tags = cmds
        .iter()
        .flat_map(|c| c.tags().iter().cloned())
        .collect::<BTreeSet<_>>();
    (keys.into_iter().collect(), tags.into_iter().collect())
}

fn header(keys: &[String], tags: &[String]) -> Vec<String> {
    let mut header: Vec<String> = vec!["step".into(), "cmd".into()];
    header.extend(keys.iter().cloned());
    header.extend(tags.iter().cloned());
    header.extend(
        [
            "wall(s)", "user(s)", "sys(s)", "max_rss", "exit", "status", "argv",
        ]
        .iter()
        .map(|s| s.to_string()),
    );
    header
}

/// The step's own line: measured wall clock, and its commands' CPU added up.
fn step_row(step: &Step, dimensions: usize) -> Vec<Cell> {
    let mut cells = vec![Cell::left(step.label()), Cell::left("-")];
    cells.extend(std::iter::repeat_n(Cell::left("-"), dimensions));

    cells.push(Cell::left(match step.wall_s() {
        Some(wall) => format!("{wall:.2}"),
        None => "-".to_string(),
    }));

    let timings: Vec<&Timing> = step
        .cmds()
        .iter()
        .filter_map(|c| match c.status() {
            // a command killed on its deadline still burned everything it says
            // it burned, so it counts toward the step
            Status::Finished(t) | Status::TimedOut(t) => Some(t),
            _ => None,
        })
        .collect();

    // no peak means nothing finished, so there is nothing to add up either
    match timings.iter().map(|t| t.max_rss_kb).max() {
        None => cells.extend(std::iter::repeat_n(Cell::left("-"), 3)),
        Some(peak) => {
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
            cells.push(Cell::left(format_bytes(peak)));
        }
    }

    // exit and status belong to commands; a count or a rollup here would be a
    // different quantity sharing a column
    cells.extend(std::iter::repeat_n(Cell::left("-"), 3));
    cells
}

/// One command's line. `first` is the step name for a collapsed step of one,
/// and a right-aligned `|` or `||` otherwise.
fn cmd_row(first: Cell, cmd: &Cmd, keys: &[String], tags: &[String]) -> Vec<Cell> {
    let mut cells = vec![first, Cell::left(cmd.label())];
    cells.extend(
        keys.iter()
            .map(|k| Cell::left(cmd.fields().get(k).map(String::as_str).unwrap_or("-"))),
    );
    cells.extend(
        tags.iter()
            .map(|t| Cell::left(if cmd.tags().contains(t) { "x" } else { "-" })),
    );

    // two separate questions: what it cost, and how it went. a command that
    // could not start has nothing to say about the first
    match cmd.status() {
        Status::Finished(t) | Status::TimedOut(t) => {
            cells.push(Cell::left(format!("{:.2}", t.wall_s)));
            cells.push(Cell::left(format!("{:.2}", t.user_s)));
            cells.push(Cell::left(format!("{:.2}", t.sys_s)));
            cells.push(Cell::left(format_bytes(t.max_rss_kb)));
            cells.push(Cell::left(t.exit.to_string()));
        }
        _ => cells.extend(std::iter::repeat_n(Cell::left("-"), 5)),
    }

    cells.push(Cell::left(status_word(cmd.status())));
    cells.push(Cell::left(cmd.line()));
    cells
}

/// The status column. It answers one question — did this work — and leaves the
/// exit column to say what the command actually reported, which is also how
/// "never started" tells itself apart from "started and failed" without a word
/// of its own: there is no exit code beside it.
fn status_word(status: &Status) -> &'static str {
    match status {
        Status::NotRun => "-",
        Status::Skipped => "skip",
        Status::TimedOut(_) => "time",
        Status::Finished(t) if t.ok() => "ok",
        _ => "fail",
    }
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

/// Lay out the rows, every column padded to its widest cell, under a commented
/// header and dashed separator if `show_header` says so.
///
/// The header sets the width of every column whether it is printed or not. That
/// is what keeps blocks lined up under a header printed once at the top: they
/// share its widths as a floor, and only drift apart where a value is wider than
/// the label above it.
fn render(header: &[String], rows: &[Vec<Cell>], show_header: bool, floor: &[usize]) -> String {
    let head = head_cells(header);

    // a floor from a different set of columns is no floor at all, which is also
    // how the modes that pass an empty one opt out
    let mut widths = if floor.len() == header.len() {
        floor.to_vec()
    } else {
        vec![0; header.len()]
    };
    widen(&mut widths, &head);
    for row in rows {
        widen(&mut widths, row);
    }

    let mut out = String::new();
    if show_header {
        let last = header.len() - 1;
        let mut sep: Vec<Cell> = widths.iter().map(|w| Cell::left("-".repeat(*w))).collect();
        sep[0] = Cell::left(format!("# {}", "-".repeat(widths[0].saturating_sub(2))));
        // argv is unpadded, so underline the label rather than the whole column
        sep[last] = Cell::left("-".repeat(header[last].chars().count()));

        write_row(&mut out, &head, &widths);
        write_row(&mut out, &sep, &widths);
    }
    for row in rows {
        write_row(&mut out, row, &widths);
    }
    out
}

/// The header as cells. The "# " marker is absorbed into the first column's
/// width, so the labels stay lined up over the data below them.
fn head_cells(header: &[String]) -> Vec<Cell> {
    let mut head: Vec<Cell> = header.iter().map(Cell::left).collect();
    head[0] = Cell::left(format!("# {}", header[0]));
    head
}

fn widen(widths: &mut [usize], cells: &[Cell]) {
    for (i, cell) in cells.iter().enumerate() {
        widths[i] = widths[i].max(cell.text.chars().count());
    }
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
