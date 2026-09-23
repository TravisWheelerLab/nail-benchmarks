//! The shape every table in this family has, without the columns that tell
//! them apart.
//!
//! A file is a format line, a preamble, one block per shard, rows sorted by
//! (query, target) within a block, and a counted trailer. That much is the
//! same whether a row carries a score per tool or a score per run, so it is
//! read here once rather than in each table's reader.
//!
//! The file is bigger than the machine, so an analysis is a pass over it
//! rather than a structure built from it. A row is a borrowed line with its
//! fields located, and what a cell means is worked out when it is asked for:
//! a summary reads the pass string and never parses a score at all.
//!
//! What is rejected here is what would make an answer wrong rather than
//! merely odd: a legend that disagrees with the runs, a row of the wrong
//! width, a block out of order, a shard twice, a count that does not match
//! the trailer. Which format line a file must open with is left to the
//! caller.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
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

pub struct Frame<R> {
    src: BufReader<R>,
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
        let mut src = BufReader::with_capacity(BUFFER, src);
        let mut buf = Vec::new();
        let mut at = 0u64;

        let mut preamble = Preamble::default();
        let mut format = String::new();
        let shard_at;

        loop {
            let Some(line) = read(&mut src, &mut buf, &mut at)? else {
                bail!("{file} ends before its header");
            };

            let text =
                std::str::from_utf8(line).with_context(|| format!("{file}:{at} is not text"))?;

            // the first line says which table this is, and nothing else may
            // come before it: a file that opens any other way was written in
            // a shape no reader here parses
            if at == 1 {
                ensure!(
                    text.starts_with("#= format "),
                    "{file} does not open a `#= format` line; \
                     it was written in an older shape",
                );

                format = text.to_string();
                continue;
            }

            let Some(rest) = text.strip_prefix("#=") else {
                ensure!(
                    text.starts_with('#'),
                    "{file}:{at} is a row before any `#= shard` marker"
                );
                continue;
            };

            let mut fields = rest.split_whitespace();
            let Some(key) = fields.next() else { continue };
            let fields: Vec<&str> = fields.collect();

            // the first block marker ends the preamble, and the shard it
            // names is the one the rows after it are in
            if key == "shard" {
                shard_at = shard(&fields)?;
                break;
            }

            preamble
                .absorb(key, &fields)
                .with_context(|| format!("{file}:{at}"))?;
        }

        let meta = preamble.finish().with_context(|| format!("in {file}"))?;

        Ok(Frame {
            src,
            buf,
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
            line: at,
            done: false,
            verdict: None,
            cutoff: (None, None),
            cutoff_for: Vec::new(),
        })
    }

    /// The `#= format` line this file opened with.
    pub fn format(&self) -> &str {
        &self.format
    }

    pub fn file(&self) -> &str {
        &self.file
    }

    /// Fix where the columns are, once the preamble has declared how many runs and
    /// tools the file carries. Rows cannot be read before this.
    pub fn layout(&mut self, layout: Layout) {
        self.layout = Some(layout);
    }

    /// Fix how a row is asked whether a run cleared its family's cutoff.
    pub fn verdict(&mut self, verdict: Verdict) {
        self.verdict = Some(verdict);
    }

    /// Advance to the next row. `false` at the trailer.
    pub fn step(&mut self) -> anyhow::Result<bool> {
        loop {
            if self.done {
                return Ok(false);
            }

            let Some(line) = read(&mut self.src, &mut self.buf, &mut self.line)? else {
                bail!(
                    "{} ends without a `#= end` line; the run that wrote it did not finish",
                    self.file
                );
            };

            if line.first() == Some(&b'#') {
                let text = std::str::from_utf8(line)
                    .with_context(|| format!("{}:{} is not text", self.file, self.line))?;

                let Some(rest) = text.strip_prefix("#=") else {
                    continue;
                };

                let mut fields = rest.split_whitespace();
                let Some(key) = fields.next() else { continue };
                let fields: Vec<&str> = fields.collect();

                match key {
                    "shard" => {
                        let shard = shard(&fields)?;
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
                        let count: u64 = fields
                            .first()
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

            if line.is_empty() {
                continue;
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

        self.spans.clear();

        let mut at = 0usize;
        while at < self.buf.len() {
            while at < self.buf.len() && self.buf[at].is_ascii_whitespace() {
                at += 1;
            }

            if at == self.buf.len() {
                break;
            }

            let start = at;
            while at < self.buf.len() && !self.buf[at].is_ascii_whitespace() {
                at += 1;
            }

            self.spans.push((start, at));
        }

        let want = layout.fields;
        ensure!(
            self.spans.len() == want,
            "{}:{} has {} fields, the header names {want}",
            self.file,
            self.line,
            self.spans.len()
        );

        if let Some(Verdict { cutoffs, .. }) = &self.verdict {
            let (start, end) = self.spans[0];

            if self.buf[start..end] != self.cutoff_for[..] {
                self.cutoff_for.clear();
                self.cutoff_for.extend_from_slice(&self.buf[start..end]);

                self.cutoff = match std::str::from_utf8(&self.buf[start..end]) {
                    Ok(family) => cutoffs.pair(family),
                    Err(_) => (None, None),
                };
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


    /// Which block the current row is in -- the unit it was searched against.
    pub fn shard(&self) -> &str {
        &self.shard
    }


    /// Whether one run reported this pair at or above its family's cutoff.
    pub fn passed(&self, run: usize) -> bool {
        let verdict = self
            .verdict
            .as_ref()
            .expect("a frame is given its verdict before its first row");

        let Some(score) = crate::scores::scan::score(self.field(verdict.at + run)) else {
            return false;
        };

        let cutoff = match self.meta.runs[run].tool {
            Tool::Nail | Tool::Hmmer => self.cutoff.0,
            Tool::Mmseqs => self.cutoff.1,
        };

        cutoff.is_some_and(|cutoff| score >= cutoff)
    }
}

/// One line, without its newline, into a buffer the reader keeps.
fn read<'a, R: Read>(
    src: &mut BufReader<R>,
    buf: &'a mut Vec<u8>,
    line: &mut u64,
) -> std::io::Result<Option<&'a [u8]>> {
    buf.clear();

    if src.read_until(b'\n', buf)? == 0 {
        return Ok(None);
    }

    *line += 1;

    while let Some(b'\n' | b'\r') = buf.last() {
        buf.pop();
    }

    Ok(Some(buf))
}

fn shard(fields: &[&str]) -> anyhow::Result<String> {
    match *fields.first().context("a `#= shard` line names no shard")? {
        "-" => Ok(String::new()),
        shard => Ok(shard.to_string()),
    }
}
