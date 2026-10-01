//! The shape every table in this family has, without the columns that tell
//! them apart.
//!
//! A file is a format line, a preamble, one block per shard, rows sorted by
//! (query, target) within a block, and a counted trailer. That much is the
//! same whether a row carries a score per tool or a score per run, so it is
//! read here once rather than in each table's reader.
//!
//! The file is bigger than the machine, so an analysis is a pass over it
//! rather than a structure built from it. toil's `Reader` cuts each line into
//! cells, and [`Frame::each`] hands a caller each row still borrowed from the
//! reader's buffer. What a cell means is worked out when it is asked for: a
//! summary reads the pass string and never parses a score at all.
//!
//! What is rejected here is what would make an answer wrong rather than
//! merely odd: a legend that disagrees with the runs, a row of the wrong
//! width, a block out of order, a shard twice, a count that does not match
//! the trailer. Which format line a file must open with is left to the
//! caller.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use anyhow::{Context, bail, ensure};

use crate::scores::{Meta, Named, Preamble, Tool};

/// How much of the file is held at once.
const BUFFER: usize = 4 << 20;

/// Where the columns every table here carries sit in a row.
pub struct Layout {
    /// How many fields a row holds.
    //
    // the score columns in the middle are the only part that
    // moves between tables; each table computes its own start,
    // and what is here is only what the envelope checks
    pub fields: usize,
}

/// Where the score columns start, and the cutoffs a score is held against.
//
// a table with a score column per run has the comparison in it already, and a
// character per run repeating the answer is a column per run of it
pub struct Verdict {
    pub at: usize,
    pub cutoffs: Named,
}

pub struct Frame<R: Read> {
    reader: toil::Reader<BufReader<R>>,
    pub meta: Meta,
    state: State,
}

/// What a frame holds rows against, and what it has seen of the file so far.
struct State {
    /// Which query and target the last row named, as `query\0target`, which
    /// orders as the pair does.
    prev: Vec<u8>,
    /// The same, being built for the row in hand.
    next: Vec<u8>,

    /// The `#= format` line the file opened with, for the caller to check.
    format: String,
    /// What to call the file in an error.
    file: String,
    layout: Option<Layout>,
    verdict: Option<Verdict>,

    /// The nail and mmseqs cutoffs of the family the row in hand names.
    //
    // a query's rows are adjacent, so this is one lookup per family
    // rather than one per cell, which at 164 runs is 164 times fewer
    cutoff: (Option<f32>, Option<f32>),
    /// Which family those are for.
    cutoff_for: Vec<u8>,

    shard: String,
    seen: HashSet<String>,
    rows: u64,
    line: u64,
    done: bool,
}

/// One row, its cells borrowed from the reader for as long as a caller looks
/// at it.
pub struct Row<'a> {
    cells: toil::Cells<'a>,
    meta: &'a Meta,
    state: &'a State,
}

impl Frame<File> {
    pub fn open(path: &Path) -> anyhow::Result<Frame<File>> {
        let file =
            File::open(path).with_context(|| format!("failed to open {}", path.display()))?;

        Frame::new(file, &path.display().to_string())
    }
}

impl<R: Read> Frame<R> {
    /// Read the format line, the preamble and the first block marker, leaving
    /// the first row next.
    pub fn new(src: R, file: &str) -> anyhow::Result<Frame<R>> {
        let mut reader = toil::Reader::new(BufReader::with_capacity(BUFFER, src))
            .with_context(|| format!("failed to read {file}"))?;

        // the first `#=` line says which table this is: a file that opens any
        // other way was written in a shape no reader here parses
        let format = reader
            .header()
            .meta_rows()
            .next()
            .filter(|row| row.key() == "format")
            .and_then(|row| row.rest(0))
            .with_context(|| {
                format!("{file} does not open a `#= format` line; it was written in an older shape")
            })?
            .to_string();

        let mut preamble = Preamble::default();
        for row in reader
            .header()
            .meta_rows()
            .filter(|row| row.key() != "format")
        {
            preamble
                .absorb(&row)
                .with_context(|| format!("in the preamble of {file}"))?;
        }

        // the first block marker ends the preamble, and the shard it names is
        // the one the rows after it are in
        let (shard_at, at) = loop {
            match reader.next_entry()? {
                None => bail!("{file} ends before its first `#= shard` marker"),
                Some((at, toil::Entry::Row(_))) => {
                    bail!("{file}:{at} is a row before any `#= shard` marker")
                }
                Some((_, toil::Entry::Comment(_))) => {}
                Some((at, toil::Entry::Meta(row))) => match row.key() {
                    "shard" => break (row.get(0).unwrap_or_default().to_string(), at),
                    _ => preamble
                        .absorb(&row)
                        .with_context(|| format!("{file}:{at}"))?,
                },
            }
        };

        let meta = preamble.finish().with_context(|| format!("in {file}"))?;

        Ok(Frame {
            reader,
            meta,
            state: State {
                prev: Vec::new(),
                next: Vec::new(),
                format,
                file: file.to_string(),
                layout: None,
                verdict: None,
                cutoff: (None, None),
                cutoff_for: Vec::new(),
                seen: HashSet::from([shard_at.clone()]),
                shard: shard_at,
                rows: 0,
                line: at as u64,
                done: false,
            },
        })
    }

    /// What the `#= format` line this file opened with says, after the key.
    pub fn format(&self) -> &str {
        &self.state.format
    }

    pub fn file(&self) -> &str {
        &self.state.file
    }

    /// Fix where the columns are, once the preamble has declared how many runs and
    /// tools the file carries. Rows cannot be read before this.
    pub fn layout(&mut self, layout: Layout) {
        self.state.layout = Some(layout);
    }

    /// Fix how a row is asked whether a run cleared its family's cutoff.
    pub fn verdict(&mut self, verdict: Verdict) {
        self.state.verdict = Some(verdict);
    }

    /// Hand every row to `f` in file order, then check the trailer.
    pub fn each(&mut self, mut f: impl FnMut(&Row) -> anyhow::Result<()>) -> anyhow::Result<()> {
        let Frame {
            reader,
            meta,
            state,
        } = self;

        while !state.done {
            let Some((at, entry)) = reader.next_entry()? else {
                bail!(
                    "{} ends without a `#= end` line; the run that wrote it did not finish",
                    state.file
                );
            };
            state.line = at as u64;

            match entry {
                toil::Entry::Comment(_) => {}
                toil::Entry::Meta(row) => state.marker(&row)?,
                toil::Entry::Row(cells) => {
                    if cells.is_empty() {
                        continue;
                    }

                    state.check(&cells)?;
                    state.rows += 1;

                    f(&Row { cells, meta, state })?;
                }
            }
        }

        Ok(())
    }
}

impl State {
    /// Take in a `#=` line between the rows: a block's marker or the trailer.
    fn marker(&mut self, row: &toil::MetaRow) -> anyhow::Result<()> {
        match row.key() {
            "shard" => {
                let shard = row.get(0).unwrap_or_default().to_string();
                ensure!(
                    self.seen.insert(shard.clone()),
                    "{}:{} opens shard {shard:?} a second time",
                    self.file,
                    self.line
                );

                self.shard = shard;
                // a target lives in one shard, so the order starts again at
                // every block
                self.prev.clear();
            }
            "end" => {
                let count: u64 = row
                    .get(0)
                    .context("a `#= end` line wants a row count")?
                    .parse()
                    .with_context(|| format!("{}:{}", self.file, self.line))?;

                ensure!(
                    count == self.rows,
                    "{} says it holds {count} rows and holds {}",
                    self.file,
                    self.rows
                );

                self.done = true;
            }
            other => bail!("{}:{} has a `#= {other}` line", self.file, self.line),
        }

        Ok(())
    }

    /// Check a row against what the file declared.
    fn check(&mut self, cells: &toil::Cells) -> anyhow::Result<()> {
        let layout = self
            .layout
            .as_ref()
            .expect("a frame is given its layout before its first row");

        let want = layout.fields;
        ensure!(
            cells.len() == want,
            "{}:{} has {} fields, the header names {want}",
            self.file,
            self.line,
            cells.len()
        );

        let (query, target) = (field(cells, 0), field(cells, 1));

        if let Some(Verdict { cutoffs, .. }) = &self.verdict
            && query != &self.cutoff_for[..]
        {
            self.cutoff_for.clear();
            self.cutoff_for.extend_from_slice(query);

            self.cutoff = match std::str::from_utf8(query) {
                Ok(family) => cutoffs.pair(family),
                Err(_) => (None, None),
            };
        }

        // the pair as `query\0target`, which orders as the pair does and costs
        // a copy of thirty bytes rather than two strings
        self.next.clear();
        self.next.extend_from_slice(query);
        self.next.push(0);
        self.next.extend_from_slice(target);

        if !self.prev.is_empty() {
            ensure!(
                self.next > self.prev,
                "{}:{} is out of order in shard {:?}",
                self.file,
                self.line,
                self.shard
            );
        }

        std::mem::swap(&mut self.prev, &mut self.next);
        Ok(())
    }
}

/// Cell `at` of a row whose width has been checked against the layout.
fn field<'a>(cells: &toil::Cells<'a>, at: usize) -> &'a [u8] {
    cells
        .field(at)
        .expect("the row's width is checked against the layout")
}

impl<'a> Row<'a> {
    /// One field of the row, by its place in the line.
    pub fn field(&self, at: usize) -> &'a [u8] {
        field(&self.cells, at)
    }

    /// Which block the row is in -- the unit it was searched against.
    pub fn shard(&self) -> &'a str {
        &self.state.shard
    }

    /// Whether one run reported this pair at or above its family's cutoff.
    pub fn passed(&self, run: usize) -> bool {
        let verdict = self
            .state
            .verdict
            .as_ref()
            .expect("a frame is given its verdict before its first row");

        let Some(score) = crate::scores::scan::score(self.field(verdict.at + run)) else {
            return false;
        };

        let cutoff = match self.meta.runs[run].tool {
            Tool::Nail => self.state.cutoff.0,
            Tool::Mmseqs => self.state.cutoff.1,
        };

        cutoff.is_some_and(|cutoff| score >= cutoff)
    }
}
