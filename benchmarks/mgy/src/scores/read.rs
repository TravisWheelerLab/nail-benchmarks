//! Reading recall's `scores.tbl` back, a row at a time.
//!
//! The file is bigger than the machine, so an analysis is a pass over it
//! rather than a structure built from it. A row is a borrowed line with its
//! fields located, and what a cell means is worked out when it is asked for:
//! a summary reads the pass string and never parses a score at all.
//!
//! What the reader refuses is what would make an answer wrong rather than
//! merely odd: another format, a legend that disagrees with the runs, a row
//! of the wrong width, a block out of order, a shard twice, a count that
//! does not match the trailer.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{Context, bail, ensure};

use super::{FORMAT, Meta, Preamble, SIGNIFICANT, Tool};

/// How much of the file is held at once.
const BUFFER: usize = 4 << 20;

pub struct Reader<R> {
    src: BufReader<R>,
    buf: Vec<u8>,
    spans: Vec<(usize, usize)>,

    /// Which query and target the last row named, as `query\0target`, which
    /// orders as the pair does.
    prev: Vec<u8>,
    /// The same, being built for the row in hand.
    next: Vec<u8>,

    pub meta: Meta,
    tools: Vec<Tool>,
    /// What to call the file in an error.
    file: String,

    shard: String,
    seen: HashSet<String>,
    rows: u64,
    line: u64,
    done: bool,
}

/// One row, as far as it has been read. A cell is parsed when it is asked for.
pub struct Row<'a> {
    line: &'a [u8],
    spans: &'a [(usize, usize)],
    tools: &'a [Tool],
}

impl Reader<File> {
    pub fn open(path: &Path) -> anyhow::Result<Reader<File>> {
        let file =
            File::open(path).with_context(|| format!("failed to open {}", path.display()))?;

        Reader::new(file, &path.display().to_string())
    }
}

impl<R: Read> Reader<R> {
    /// Read the preamble and the header, leaving the first row next.
    pub fn new(src: R, file: &str) -> anyhow::Result<Reader<R>> {
        let mut src = BufReader::with_capacity(BUFFER, src);
        let mut buf = Vec::new();
        let mut at = 0u64;

        let mut preamble = Preamble::default();
        let shard_at;

        loop {
            let Some(line) = read(&mut src, &mut buf, &mut at)? else {
                bail!("{file} ends before its header");
            };

            if at == 1 {
                ensure!(
                    line == FORMAT.as_bytes(),
                    "{file} does not open `{FORMAT}`; it was written in an older shape",
                );
                continue;
            }

            let text = std::str::from_utf8(line)
                .with_context(|| format!("{file}:{at} is not text"))?;

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
        let tools = meta.tools();

        Ok(Reader {
            src,
            buf,
            spans: Vec::new(),
            prev: Vec::new(),
            next: Vec::new(),
            meta,
            tools,
            file: file.to_string(),
            seen: HashSet::from([shard_at.clone()]),
            shard: shard_at,
            rows: 0,
            line: at,
            done: false,
        })
    }

    /// The next row, or `None` at the end of the file.
    pub fn next(&mut self) -> anyhow::Result<Option<Row<'_>>> {
        loop {
            if self.done {
                return Ok(None);
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
                        return Ok(None);
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

            return Ok(Some(Row {
                line: &self.buf,
                spans: &self.spans,
                tools: &self.tools,
            }));
        }
    }

    /// Locate a row's fields and check it against what the file declared.
    fn row(&mut self) -> anyhow::Result<()> {
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

        // query, target, pass, one per tool, inc, dom
        let want = 5 + self.tools.len();
        ensure!(
            self.spans.len() == want,
            "{}:{} has {} fields, the header names {want}",
            self.file,
            self.line,
            self.spans.len()
        );

        let pass = self.field(2);

        ensure!(
            pass.len() == self.meta.runs.len(),
            "{}:{} has a pass string of {} for {} runs",
            self.file,
            self.line,
            pass.len(),
            self.meta.runs.len()
        );

        for (at, letter) in pass.iter().enumerate() {
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

    fn field(&self, at: usize) -> &[u8] {
        let (start, end) = self.spans[at];
        &self.buf[start..end]
    }

}

impl Row<'_> {
    fn field(&self, at: usize) -> &[u8] {
        let (start, end) = self.spans[at];
        &self.line[start..end]
    }

    pub fn query(&self) -> &str {
        text(self.field(0))
    }

    pub fn target(&self) -> &str {
        text(self.field(1))
    }

    /// One character per run, in the order `#= pass` names them.
    pub fn pass(&self) -> &[u8] {
        self.field(2)
    }

    /// Whether one run reported this pair at or above its family's cutoff.
    pub fn passed(&self, run: usize) -> bool {
        self.pass()
            .get(run)
            .is_some_and(u8::is_ascii_uppercase)
    }

    /// Whether any run of this tool reported the pair, at any score.
    pub fn present(&self, tool: Tool) -> bool {
        self.at(tool).is_some_and(|at| self.field(at) != b"-")
    }

    /// What this tool scored the pair, whichever of its runs reported it.
    pub fn score(&self, tool: Tool) -> Option<f32> {
        let at = self.at(tool)?;
        super::scan::score(self.field(at))
    }

    /// hmmer's inclusion count.
    pub fn inc(&self) -> Option<u32> {
        let at = 3 + self.tools.len();
        text(self.field(at)).parse().ok()
    }

    /// hmmer's domain scores, in the order the domtbl listed them.
    pub fn domains(&self) -> impl Iterator<Item = f32> + '_ {
        let at = 4 + self.tools.len();
        let field = self.field(at);

        field
            .split(|byte| *byte == b',')
            .filter_map(|score| super::scan::score(score))
    }

    /// How many of the pair's domains carry enough of the hit to count as
    /// their own.
    ///
    /// The measure is against the best domain rather than an absolute score:
    /// what is being asked is whether the hit is one region of the sequence or
    /// several, and a weak family's several are still several.
    pub fn domain_count(&self) -> usize {
        let best = self.domains().fold(f32::NEG_INFINITY, f32::max);
        if best == f32::NEG_INFINITY {
            return 0;
        }

        match best > 0.0 {
            true => self
                .domains()
                .filter(|score| score / best >= SIGNIFICANT)
                .count(),
            false => self.domains().count(),
        }
    }

    /// Where a tool's score column sits in this row, if the table carries one.
    fn at(&self, tool: Tool) -> Option<usize> {
        self.tools
            .iter()
            .position(|carried| *carried == tool)
            .map(|at| 3 + at)
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

/// A cell as text. Every cell this table writes is ascii, and a reader that
/// got here has already had the line checked.
fn text(field: &[u8]) -> &str {
    std::str::from_utf8(field).unwrap_or("")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two runs, two rows, and a pair hmmer never reported.
    const FILE: &str = "\
#= format scores 2
#= query 3 18 120
#= target 1 3 12 40
#= cutoffs /x/cutoffs.tbl c=0
#= run nail-a nail 3.0000 s=9.0
#= run hmmer hmmer 2.0000
#= pass nail-a hmmer
# query target           pass nail   hmmer  inc dom
# ----- ---------------- ---- ------ ------ --- ---
#= shard 1
alpha MGYP000000000001 NH 25.0   24.0   1   -3.0,2.0
beta MGYP000000000002 Nh 11.0   -      -   -
#= end 2
";

    fn read(text: &str) -> anyhow::Result<Vec<String>> {
        let mut reader = Reader::new(text.as_bytes(), "test")?;
        let mut out = Vec::new();

        while let Some(row) = reader.next()? {
            out.push(format!(
                "{} {} {} {:?} {:?} {:?} {} {}",
                row.query(),
                row.target(),
                String::from_utf8_lossy(row.pass()),
                row.score(Tool::Nail),
                row.score(Tool::Hmmer),
                row.inc(),
                row.present(Tool::Hmmer),
                row.domain_count(),
            ));
        }

        Ok(out)
    }

    #[test]
    fn a_row_is_read_as_it_was_written() {
        let rows = read(FILE).unwrap();

        assert_eq!(
            rows,
            [
                "alpha MGYP000000000001 NH Some(25.0) Some(24.0) Some(1) true 1",
                "beta MGYP000000000002 Nh Some(11.0) None None false 0",
            ]
        );
    }

    #[test]
    fn the_domain_count_does_not_depend_on_the_order() {
        let swapped = FILE.replace("-3.0,2.0", "2.0,-3.0");
        assert_eq!(read(&swapped).unwrap(), read(FILE).unwrap());
    }

    #[test]
    fn a_run_is_read_by_where_it_sits() {
        let mut reader = Reader::new(FILE.as_bytes(), "test").unwrap();
        let row = reader.next().unwrap().unwrap();

        assert!(row.passed(0));
        assert!(row.passed(1));

        let row = reader.next().unwrap().unwrap();
        assert!(row.passed(0));
        assert!(!row.passed(1));
    }

    /// Every way the file can be wrong that would otherwise be read as an
    /// answer rather than as a fault.
    #[test]
    fn a_broken_file_is_refused() {
        let cases = [
            ("an older format", FILE.replace("#= format scores 2\n", "")),
            (
                "a legend that names another run",
                FILE.replace("#= pass nail-a hmmer", "#= pass nail-a mm"),
            ),
            (
                "a pass string of the wrong length",
                FILE.replace(" NH ", " N "),
            ),
            (
                "a letter that is not the run's tool",
                FILE.replace(" NH ", " NM "),
            ),
            (
                "a row of the wrong width",
                FILE.replace("11.0   -      -   -", "11.0   -      -"),
            ),
            (
                "a shard opened twice",
                FILE.replace("#= end 2", "#= shard 1\n#= end 2"),
            ),
            (
                "rows out of order",
                FILE.replace(
                    "alpha MGYP000000000001 NH 25.0   24.0   1   -3.0,2.0\nbeta MGYP000000000002 Nh 11.0   -      -   -",
                    "beta MGYP000000000002 Nh 11.0   -      -   -\nalpha MGYP000000000001 NH 25.0   24.0   1   -3.0,2.0",
                ),
            ),
            ("no trailer", FILE.replace("#= end 2\n", "")),
            ("a trailer that counts wrong", FILE.replace("#= end 2", "#= end 3")),
            (
                "a row before any shard",
                FILE.replace("#= shard 1\n", ""),
            ),
        ];

        for (what, text) in cases {
            assert!(read(&text).is_err(), "{what} was read without complaint");
        }
    }
}
