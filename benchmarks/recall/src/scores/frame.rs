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
//! cells, and a row is those cells copied into one buffer reused for every
//! row. What a cell means is worked out when it is asked for: a summary reads
//! the pass string and never parses a score at all.
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

use crate::scores::{Meta, Preamble, SIGNIFICANT, Tool};

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
    /// Where the pass string sits.
    pub pass: usize,
    pub dom: usize,
}

pub struct Frame<R: Read> {
    reader: toil::Reader<BufReader<R>>,
    /// The row in hand, copied out of the reader so a caller can ask for its
    /// fields after `step` returns.
    buf: Vec<u8>,
    spans: Vec<(usize, usize)>,

    /// Which query and target the last row named, as `query\0target`, which
    /// orders as the pair does.
    prev: Vec<u8>,
    /// The same, being built for the row in hand.
    next: Vec<u8>,

    pub meta: Meta,
    /// The `#= format` line the file opened with, for the caller to check.
    format: String,
    /// What to call the file in an error.
    file: String,
    layout: Option<Layout>,

    shard: String,
    seen: HashSet<String>,
    rows: u64,
    line: u64,
    done: bool,
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

        // the first line says which table this is, and nothing else may come
        // before it: a file that opens any other way was written in a shape no
        // reader here parses
        let format = reader
            .header()
            .preamble()
            .first()
            .filter(|line| line.starts_with("#= format "))
            .with_context(|| {
                format!("{file} does not open a `#= format` line; it was written in an older shape")
            })?
            .clone();

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
            buf: Vec::new(),
            spans: Vec::new(),
            prev: Vec::new(),
            next: Vec::new(),
            meta,
            format,
            file: file.to_string(),
            layout: None,
            seen: HashSet::from([shard_at.clone()]),
            shard: shard_at,
            rows: 0,
            line: at as u64,
            done: false,
        })
    }

    /// The `#= format` line this file opened with.
    pub fn format(&self) -> &str {
        &self.format
    }

    /// Fix where the columns are, once the preamble has declared how many runs and
    /// tools the file carries. Rows cannot be read before this.
    pub fn layout(&mut self, layout: Layout) {
        self.layout = Some(layout);
    }

    /// Advance to the next row. `false` at the trailer.
    pub fn step(&mut self) -> anyhow::Result<bool> {
        loop {
            if self.done {
                return Ok(false);
            }

            let Some((at, entry)) = self.reader.next_entry()? else {
                bail!(
                    "{} ends without a `#= end` line; the run that wrote it did not finish",
                    self.file
                );
            };
            self.line = at as u64;

            let cells = match entry {
                toil::Entry::Row(cells) => cells,
                toil::Entry::Comment(_) => continue,
                toil::Entry::Meta(row) => {
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
                            // a target lives in one shard, so the order starts
                            // again at every block
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
                            return Ok(false);
                        }
                        other => bail!("{}:{} has a `#= {other}` line", self.file, self.line),
                    }
                    continue;
                }
            };
            if cells.is_empty() {
                continue;
            }

            self.buf.clear();
            self.spans.clear();
            for i in 0..cells.len() {
                let start = self.buf.len();
                self.buf
                    .extend_from_slice(cells.field(i).expect("i is below the cell count"));
                self.spans.push((start, self.buf.len()));
            }

            self.row()?;
            self.rows += 1;

            return Ok(true);
        }
    }

    /// Locate a row's fields and check it against what the file declared.
    fn row(&mut self) -> anyhow::Result<()> {
        let layout = self
            .layout
            .as_ref()
            .expect("a frame is given its layout before its first row");

        let want = layout.fields;
        ensure!(
            self.spans.len() == want,
            "{}:{} has {} fields, the header names {want}",
            self.file,
            self.line,
            self.spans.len()
        );

        {
            let (start, end) = self.spans[layout.pass];
            let pass = &self.buf[start..end];

            ensure!(
                pass.len() == self.meta.runs.len(),
                "{}:{} has a pass string of {} for {} runs",
                self.file,
                self.line,
                pass.len(),
                self.meta.runs.len()
            );

            for (at, letter) in pass.iter().enumerate() {
                // `-` is the one character that is not a tool's: it says the
                // run's seeding never offered the pair, which is a different
                // answer from the run having scored it badly
                if *letter == b'-' {
                    continue;
                }

                let tool = Tool::of_letter(*letter);
                ensure!(
                    tool == Some(self.meta.runs[at].tool),
                    "{}:{} has {:?} where run {} is {}",
                    self.file,
                    self.line,
                    *letter as char,
                    self.meta.runs[at].name,
                    self.meta.runs[at].tool
                );
            }
        }

        // the pair as `query\0target`, which orders as the pair does and costs
        // a copy of thirty bytes rather than two strings
        let (query, target) = (self.spans[0], self.spans[1]);

        self.next.clear();
        self.next.extend_from_slice(&self.buf[query.0..query.1]);
        self.next.push(0);
        self.next.extend_from_slice(&self.buf[target.0..target.1]);

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

    /// One field of the row in hand, by its place in the line.
    pub fn field(&self, at: usize) -> &[u8] {
        let (start, end) = self.spans[at];
        &self.buf[start..end]
    }

    /// One character per run, in the order `#= pass` names them.
    pub fn pass(&self) -> &[u8] {
        self.field(self.at().pass)
    }

    /// Whether one run reported this pair at or above its family's cutoff.
    pub fn passed(&self, run: usize) -> bool {
        self.pass().get(run).is_some_and(u8::is_ascii_uppercase)
    }

    /// hmmer's domain scores, in the order the domtbl listed them.
    pub fn domains(&self) -> impl Iterator<Item = f32> + '_ {
        self.field(self.at().dom)
            .split(|byte| *byte == b',')
            .filter_map(super::scan::score)
    }

    /// How many of the pair's domains carry enough of the hit to count as
    /// their own.
    pub fn domain_count(&self) -> usize {
        let best = self.domains().fold(f32::NEG_INFINITY, f32::max);
        if best == f32::NEG_INFINITY {
            return 0;
        }

        // against the best domain rather than an absolute score:
        // the question is whether the hit is one region of the
        // sequence or several, and a weak family's several are
        // still several
        match best > 0.0 {
            true => self
                .domains()
                .filter(|score| score / best >= SIGNIFICANT)
                .count(),
            false => self.domains().count(),
        }
    }

    fn at(&self) -> &Layout {
        self.layout
            .as_ref()
            .expect("a frame is given its layout before its first row")
    }
}
