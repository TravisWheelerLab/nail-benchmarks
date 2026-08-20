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
use crate::fmt::{bytes, cpu_pct, dash, secs};
use crate::sink::Sink;
use crate::step::{Step, Strategy};

/// Marks a command that ran after the one above it.
const SERIAL: &str = "|";
/// Marks a command that ran alongside the others in its step.
const BATCH: &str = "||";

/// The columns every row ends with, after whatever fields and tags the run
/// carries.
const METRICS: [&str; 8] = [
    "wall(s)", "user(s)", "sys(s)", "cpu(%)", "max_rss", "exit", "status", "argv",
];

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
    /// The columns every block carries, worked out before anything runs so blocks
    /// agree whatever order fields turn up in. Unused by [`Mode::Ragged`], which
    /// asks each step instead.
    columns: Columns,
    /// Every row so far, for [`Mode::Whole`].
    rows: Vec<Vec<Cell>>,
    /// Widths every block starts from, worked out before the run from everything
    /// already known: names, fields, tags, argv. Only the numbers are missing,
    /// and their headings are wider than they usually are. `None` for the modes
    /// that do not share widths between blocks.
    floor: Option<Vec<usize>>,
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
            columns: Columns::default(),
            rows: Vec::new(),
            floor: None,
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
        self.columns = Columns::of(&cmds);

        // Ragged blocks are meant to differ, and Whole renders in one go, so
        // neither has anything to share.
        self.floor = match self.mode {
            Mode::Ragged | Mode::Whole => None,
            _ => Some(self.columns.measure(steps)),
        };

        self.text.clear();
        if matches!(
            self.mode,
            Mode::Blocks {
                headers: Headers::Once
            }
        ) {
            self.text = render(&self.columns.header(), &[], true, self.floor.as_deref());
        }
        self.flush()
    }

    // no `record`: a step's rows are built from the step itself once it is done,
    // which keeps them in the order the commands were declared rather than the
    // order a batch happened to finish them in

    fn step_done(&mut self, step: &Step) -> anyhow::Result<()> {
        let columns = match self.mode {
            Mode::Ragged => Columns::of(&step.cmds().iter().collect::<Vec<_>>()),
            _ => self.columns.clone(),
        };

        let mut rows = columns.block(step);

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
            &columns.header(),
            &rows,
            show_header,
            self.floor.as_deref(),
        ));
        self.flush()
    }

    fn finish(&mut self) -> anyhow::Result<()> {
        if self.mode == Mode::Whole {
            let rows = std::mem::take(&mut self.rows);
            self.text = render(&self.columns.header(), &rows, true, None);
            self.flush()?;
        }
        Ok(())
    }
}

/// The field keys and tags a block carries, in front of [`METRICS`].
///
/// Which ones there are depends on the commands, so every part of a block —
/// the header, the step line, each command line — has to agree about them. They
/// live here rather than being handed to each in turn.
#[derive(Clone, Debug, Default)]
struct Columns {
    keys: Vec<String>,
    tags: Vec<String>,
}

impl Columns {
    /// The keys and tags these commands carry, each in sorted order.
    fn of(cmds: &[&Cmd]) -> Columns {
        Columns {
            keys: cmds
                .iter()
                .flat_map(|c| c.fields().keys().cloned())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect(),
            tags: cmds
                .iter()
                .flat_map(|c| c.tags().iter().cloned())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect(),
        }
    }

    fn header(&self) -> Vec<String> {
        let mut header: Vec<String> = vec!["step".into(), "cmd".into()];
        header.extend(self.keys.iter().cloned());
        header.extend(self.tags.iter().cloned());
        header.extend(METRICS.iter().map(|s| s.to_string()));
        header
    }

    /// One step's rows: its own line, then a line per command.
    fn block(&self, step: &Step) -> Vec<Vec<Cell>> {
        let mut rows = Vec::new();

        // a step of one would just repeat itself, so it gets no line of its own
        // and keeps the first column instead, with its command filling in the rest
        let first = if step.cmds().len() == 1 {
            Cell::left(step.label())
        } else {
            rows.push(self.step_row(step));
            Cell::right(match step.strategy() {
                Strategy::Serial => SERIAL,
                Strategy::Batched { .. } => BATCH,
            })
        };

        for cmd in step.cmds() {
            rows.push(self.cmd_row(first.clone(), cmd));
        }
        rows
    }

    /// How wide each column has to be for every block to fit under one header.
    ///
    /// Everything but the numbers is already known before the run, and the
    /// headings above the numbers are wider than the numbers usually are.
    fn measure(&self, steps: &[Step]) -> Vec<usize> {
        let header = self.header();
        let mut widths = vec![0; header.len()];

        widen(&mut widths, &head_cells(&header));
        for step in steps {
            for row in self.block(step) {
                widen(&mut widths, &row);
            }
        }
        widths
    }
}

/// What one row has to say about what something cost and how it went. Anything
/// left out prints as `-`, which is how a command that never started says it has
/// no numbers.
#[derive(Default)]
struct Metrics {
    wall_s: Option<f64>,
    user_s: Option<f64>,
    sys_s: Option<f64>,
    max_rss_kb: Option<i64>,
    exit: Option<i32>,
    status: Option<&'static str>,
    argv: Option<String>,
}

impl Metrics {
    /// The cells, in the order [`METRICS`] names them. A column added to one
    /// without the other does not compile.
    fn cells(self) -> [Cell; METRICS.len()] {
        let cpu = match (self.user_s, self.sys_s) {
            (Some(user), Some(sys)) => cpu_pct(user + sys, self.wall_s),
            _ => dash(),
        };

        [
            Cell::left(secs(self.wall_s)),
            Cell::left(secs(self.user_s)),
            Cell::left(secs(self.sys_s)),
            Cell::left(cpu),
            Cell::left(self.max_rss_kb.map(bytes).unwrap_or_else(dash)),
            Cell::left(self.exit.map(|e| e.to_string()).unwrap_or_else(dash)),
            Cell::left(self.status.unwrap_or("-")),
            Cell::left(self.argv.unwrap_or_else(dash)),
        ]
    }
}

impl Columns {
    /// The step's own line: measured wall clock, and its commands' CPU added up.
    fn step_row(&self, step: &Step) -> Vec<Cell> {
        let mut cells = vec![Cell::left(step.label()), Cell::left("-")];
        cells.extend(std::iter::repeat_n(
            Cell::left("-"),
            self.keys.len() + self.tags.len(),
        ));

        let timings: Vec<&Timing> = step
            .cmds()
            .iter()
            .filter_map(|c| c.status().timing())
            .collect();

        // exit, status and argv belong to commands; a count or a rollup here
        // would be a different quantity sharing a column
        let mut metrics = Metrics {
            wall_s: step.wall_s(),
            ..Metrics::default()
        };

        // no peak means nothing finished, so there is nothing to add up either
        if let Some(peak) = timings.iter().map(|t| t.max_rss_kb).max() {
            // summed CPU against measured wall is what shows whether a batch
            // actually bought anything: a step that ran four at once reads about
            // four times what any one of them did
            metrics.user_s = Some(timings.iter().map(|t| t.user_s).sum());
            metrics.sys_s = Some(timings.iter().map(|t| t.sys_s).sum());
            // the largest any one process got, which is not the same as the most
            // the step held at once — wait4 cannot tell us that
            metrics.max_rss_kb = Some(peak);
        }

        cells.extend(metrics.cells());
        cells
    }

    /// One command's line. `first` is the step name for a collapsed step of one,
    /// and a right-aligned `|` or `||` otherwise.
    fn cmd_row(&self, first: Cell, cmd: &Cmd) -> Vec<Cell> {
        let mut cells = vec![first, Cell::left(cmd.label())];
        cells.extend(
            self.keys
                .iter()
                .map(|k| Cell::left(cmd.fields().get(k).map(String::as_str).unwrap_or("-"))),
        );
        cells.extend(
            self.tags
                .iter()
                .map(|t| Cell::left(if cmd.tags().contains(t) { "x" } else { "-" })),
        );

        // two separate questions: what it cost, and how it went. a command that
        // could not start has nothing to say about the first
        let t = cmd.status().timing();
        cells.extend(
            Metrics {
                wall_s: t.map(|t| t.wall_s),
                user_s: t.map(|t| t.user_s),
                sys_s: t.map(|t| t.sys_s),
                max_rss_kb: t.map(|t| t.max_rss_kb),
                exit: t.map(|t| t.exit),
                status: Some(status_word(cmd.status())),
                argv: Some(cmd.line()),
            }
            .cells(),
        );
        cells
    }
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
fn render(
    header: &[String],
    rows: &[Vec<Cell>],
    show_header: bool,
    floor: Option<&[usize]>,
) -> String {
    let head = head_cells(header);

    // a floor from a different set of columns is no floor at all
    let mut widths = match floor.filter(|floor| floor.len() == header.len()) {
        Some(floor) => floor.to_vec(),
        None => vec![0; header.len()],
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

