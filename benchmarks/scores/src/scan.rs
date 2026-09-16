//! Reading a results table without allocating per row.
//!
//! A shard's tables are a few gigabytes of short lines, read once and never
//! kept, so the cost that matters is per byte rather than per file. [`Lines`]
//! hands out a line at a time as a slice of its own buffer, and [`fields`]
//! walks a line once and stops at the last field wanted.
//!
//! Nothing here knows a layout. Which field is the query and which the score
//! belongs to `libsail`, and is read off a [`HitColumns`] at the call.

use std::io::Read;

use libsail::tbl::HitColumns;
use libsail::tbl::hmmer::HmmerTable;

/// How much of a file is held at once. Big enough that a read is sequential
/// and the readahead has something to work with.
const BUFFER: usize = 8 << 20;

/// Where a line stops being a line. A hit table's rows are tens of bytes; this
/// is here so a binary file handed to the wrong reader fails rather than
/// growing the buffer until the machine gives out.
const LONGEST: usize = 64 << 20;

/// The lines of a file that have something on them.
///
/// Blank lines and `#` lines are passed over wherever they fall, not only at
/// the top: a run split across query sets is a `cat` of its parts, so its
/// header turns up again in the middle.
pub struct Lines<R> {
    src: R,
    buf: Vec<u8>,
    /// Where the next line starts.
    at: usize,
    /// How much of the buffer holds file.
    end: usize,
    eof: bool,
    line: u64,
}

impl<R: Read> Lines<R> {
    pub fn new(src: R) -> Lines<R> {
        Lines {
            src,
            buf: vec![0; BUFFER],
            at: 0,
            end: 0,
            eof: false,
            line: 0,
        }
    }

    /// The next line with something on it, without its newline, and which
    /// line of the file it was -- the ones passed over counted, so an error
    /// can say where it was.
    pub fn next(&mut self) -> std::io::Result<Option<(u64, &[u8])>> {
        loop {
            let Some((at, end)) = self.raw()? else {
                return Ok(None);
            };

            if !blank(&self.buf[at..end]) {
                return Ok(Some((self.line, &self.buf[at..end])));
            }
        }
    }

    /// The next line, whatever is on it.
    fn raw(&mut self) -> std::io::Result<Option<(usize, usize)>> {
        loop {
            if let Some(i) = memchr::memchr(b'\n', &self.buf[self.at..self.end]) {
                let at = self.at;
                self.at = at + i + 1;
                self.line += 1;
                return Ok(Some((at, trim(&self.buf[at..at + i]) + at)));
            }

            if self.eof {
                if self.at == self.end {
                    return Ok(None);
                }

                // a last line with no newline under it
                let (at, end) = (self.at, self.end);
                self.at = self.end;
                self.line += 1;
                return Ok(Some((at, trim(&self.buf[at..end]) + at)));
            }

            self.fill()?;
        }
    }

    fn fill(&mut self) -> std::io::Result<()> {
        // what is left is the front of a line, so it moves to the front of the
        // buffer and the read goes on from there
        self.buf.copy_within(self.at..self.end, 0);
        self.end -= self.at;
        self.at = 0;

        if self.end == self.buf.len() {
            if self.buf.len() >= LONGEST {
                return Err(std::io::Error::other(format!(
                    "a line longer than {LONGEST} bytes: is this a hit table?"
                )));
            }

            self.buf.resize((self.buf.len() * 2).min(LONGEST), 0);
        }

        let read = self.src.read(&mut self.buf[self.end..])?;
        self.end += read;
        self.eof = read == 0;

        Ok(())
    }
}

/// Whether a line is one to pass over.
fn blank(line: &[u8]) -> bool {
    match line.first() {
        None | Some(b'#') => true,
        Some(_) => line.iter().all(u8::is_ascii_whitespace),
    }
}

/// How much of a line is left once a `\r` is dropped off the end.
fn trim(line: &[u8]) -> usize {
    match line.last() {
        Some(b'\r') => line.len() - 1,
        _ => line.len(),
    }
}

/// The fields at `want`, in one pass over the line.
///
/// Stops at the last index wanted, so a layout carrying what it needs in front
/// costs nothing for the rest of the row. `None` where the line has fewer
/// fields than that.
pub fn fields<const K: usize>(line: &[u8], want: [usize; K]) -> Option<[&[u8]; K]> {
    let last = *want.iter().max()?;

    let mut out = [&line[..0]; K];
    let mut got = 0usize;

    let mut field = 0usize;
    let mut at = 0usize;

    while at < line.len() {
        while at < line.len() && line[at].is_ascii_whitespace() {
            at += 1;
        }

        if at == line.len() {
            break;
        }

        let start = at;
        while at < line.len() && !line[at].is_ascii_whitespace() {
            at += 1;
        }

        for (k, &index) in want.iter().enumerate() {
            if index == field {
                out[k] = &line[start..at];
                got += 1;
            }
        }

        if field == last {
            return (got == K).then_some(out);
        }

        field += 1;
    }

    None
}

/// How many fields a line has, which is what a layout is checked against.
pub fn count(line: &[u8]) -> usize {
    line.split(|b| b.is_ascii_whitespace())
        .filter(|field| !field.is_empty())
        .count()
}

/// The query, target and score of one row, at layout `C`'s own indices.
pub fn hit<C: HitColumns>(line: &[u8]) -> Option<[&[u8]; 3]> {
    fields(line, [C::QUERY, C::TARGET, C::SCORE])
}

/// hmmer's `--tblout`, for the one column `libsail`'s layout does not carry:
/// how many of a pair's domains fall inside the inclusion threshold.
///
/// ```text
/// |  0   | 1 |   2  | 3 |  4   |  5  |  6 | ... | 15 | 16 | 17 |   18
/// # target acc query acc e-value score bias ... dom  rep  inc  description
/// ```
const INC: usize = 17;

/// The query, target, score and `inc` of one `--tblout` row.
pub fn hmmer_hit(line: &[u8]) -> Option<[&[u8]; 4]> {
    fields(
        line,
        [
            HmmerTable::QUERY,
            HmmerTable::TARGET,
            HmmerTable::SCORE,
            INC,
        ],
    )
}

/// Whether a row carries the fields layout `C` calls for.
///
/// A table's rows are all one shape, so this is asked of the first row of a
/// file and no others. What it catches is a table of the wrong tool, where
/// reading field 6 as a score would otherwise give a number.
pub fn fits<C: HitColumns>(line: &[u8]) -> bool {
    count(line) >= C::N_FIELDS
}

/// One score field, or `None` where it is not a number a cutoff can be applied
/// to -- a NaN and an infinity included, since neither compares against one in
/// a way that means anything.
pub fn score(field: &[u8]) -> Option<f32> {
    let n: f32 = std::str::from_utf8(field).ok()?.parse().ok()?;
    n.is_finite().then_some(n)
}

/// The number in an MGnify protein name: `MGYP` and twelve digits.
///
/// Twelve digits is under 2^40, and the numbers order as the names do, so a
/// target turned into one of these sorts a shard without its name being
/// interned, hashed or compared.
pub fn mgyp(name: &[u8]) -> Option<u64> {
    let digits = name.strip_prefix(b"MGYP")?;
    if digits.len() != 12 {
        return None;
    }

    let mut out = 0u64;
    for &byte in digits {
        let digit = byte.wrapping_sub(b'0');
        if digit > 9 {
            return None;
        }

        out = out * 10 + digit as u64;
    }

    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every line of a file, the plain way, to hold `Lines` against.
    fn naive(text: &str) -> Vec<String> {
        text.split('\n')
            .map(|line| line.strip_suffix('\r').unwrap_or(line))
            .filter(|line| !line.is_empty() && !line.starts_with('#') && !line.trim().is_empty())
            .map(str::to_string)
            .collect()
    }

    fn read(text: &str, buffer: usize) -> Vec<String> {
        let mut lines = Lines {
            src: text.as_bytes(),
            buf: vec![0; buffer],
            at: 0,
            end: 0,
            eof: false,
            line: 0,
        };

        let mut out = Vec::new();
        while let Some((_, line)) = lines.next().unwrap() {
            out.push(String::from_utf8(line.to_vec()).unwrap());
        }

        out
    }

    /// A buffer that ends mid-line, mid-newline, and between the two, at every
    /// size small enough to land there.
    #[test]
    fn a_line_survives_the_end_of_the_buffer() {
        let text = "a b c\n#head\n\nd e f\r\nlonger than the buffer will be\ng h i";

        for buffer in 1..40 {
            assert_eq!(read(text, buffer), naive(text), "buffer of {buffer}");
        }
    }

    #[test]
    fn a_file_reads_the_same_whole() {
        let text = "# one\nq1 t1 12.5\n\nq2 t2 -3\n# two\nq3 t3 4.0\n";
        assert_eq!(read(text, 8 << 20), naive(text));
        assert_eq!(read("", 64), Vec::<String>::new());
        assert_eq!(read("\n\n\n", 64), Vec::<String>::new());
        assert_eq!(read("q t 1", 64), vec!["q t 1"]);
    }

    #[test]
    fn a_line_too_long_is_refused() {
        let long = "x".repeat(LONGEST + 1);
        let mut lines = Lines {
            src: long.as_bytes(),
            buf: vec![0; LONGEST],
            at: 0,
            end: 0,
            eof: false,
            line: 0,
        };

        assert!(lines.next().is_err());
    }

    #[test]
    fn fields_are_what_splitting_gives() {
        let rows = [
            "7tm_1 MGYP000522683479 98.7 1.2e-30 1",
            "  leading   and    doubled   spaces  ",
            "one\ttab\tseparated\trow\there",
        ];

        for row in rows {
            let split: Vec<&str> = row.split_whitespace().collect();

            for i in 0..split.len() {
                let [got] = fields(row.as_bytes(), [i]).unwrap();
                assert_eq!(got, split[i].as_bytes(), "{row:?} field {i}");
            }

            assert_eq!(count(row.as_bytes()), split.len());
            assert!(fields(row.as_bytes(), [split.len()]).is_none());
        }
    }

    #[test]
    fn fields_come_back_in_the_order_asked_for() {
        let row = b"q t 98.7 x y 1";

        let [target, score, query] = fields(row, [1, 2, 0]).unwrap();
        assert_eq!((query, target, score), (&b"q"[..], &b"t"[..], &b"98.7"[..]));

        // the same field twice, which a layout that names one column for two
        // things would ask for
        let [a, b] = fields(row, [2, 2]).unwrap();
        assert_eq!((a, b), (&b"98.7"[..], &b"98.7"[..]));
    }

    /// A row of each of the four layouts, as the tool itself wrote it, so the
    /// indices `libsail` declares are held against a real file.
    #[test]
    fn a_hit_is_read_at_its_layout_s_indices() {
        use libsail::tbl::blast::BlastTable;
        use libsail::tbl::hmmer::HmmerDomTable;
        use libsail::tbl::nail::NailTable;

        let nail = b"MGYP005808827855 AAA_30 42     152    22    121   18.1  0.0  9.6e-5 0.066";
        assert!(fits::<NailTable>(nail));
        assert_eq!(
            hit::<NailTable>(nail).unwrap(),
            [&b"AAA_30"[..], b"MGYP005808827855", b"18.1"]
        );

        let mmseqs = b"AAA_30 MGYP005808827855 31.0 100 60 2 3 102 40 139 3.2e-12 45.7";
        assert!(fits::<BlastTable>(mmseqs));
        assert_eq!(
            hit::<BlastTable>(mmseqs).unwrap(),
            [&b"AAA_30"[..], b"MGYP005808827855", b"45.7"]
        );

        // the description field carries spaces, so the row splits into more
        // than the eighteen a layout asks for
        let hmmer = b"MGYP003392859013     -          2_5_RNA_ligase2      PF13563.10     1e-15   54.0   0.2   1.3e-15   53.7   0.2   1.2   1   0   0   1   1   1   1 CR=1 FL=0";
        assert!(fits::<HmmerTable>(hmmer));
        assert_eq!(
            hmmer_hit(hmmer).unwrap(),
            [&b"2_5_RNA_ligase2"[..], b"MGYP003392859013", b"54.0", b"1"]
        );

        let dom = b"MGYP000987338150     -            178 2-oxoacid_dh         PF00198.27   232   9.3e-44  145.5   0.5   1   1   1.1e-46   1.1e-43  145.2   0.5     5   143    42   178    39   178 0.97 FL=0";
        assert!(fits::<HmmerDomTable>(dom));
        assert_eq!(
            hit::<HmmerDomTable>(dom).unwrap(),
            [&b"2-oxoacid_dh"[..], b"MGYP000987338150", b"145.2"]
        );

        // a table of the wrong tool, which is what the field count is for
        assert!(!fits::<HmmerTable>(nail));
        assert!(!fits::<BlastTable>(hmmer_short()));
    }

    /// An hmmer row cut to the fields before the description, which is fewer
    /// than blast's twelve are wide.
    fn hmmer_short() -> &'static [u8] {
        b"MGYP003392859013 - 2_5_RNA_ligase2 PF13563.10 1e-15 54.0"
    }

    #[test]
    fn a_score_is_what_parse_gives() {
        for text in ["98.7", "-3", "0", "1e-30", "-0.0", ".5", "1.", "6.02e23"] {
            assert_eq!(score(text.as_bytes()), text.parse::<f32>().ok());
        }

        for text in ["", "-", "nan", "inf", "-inf", "98.7x", "MGYP1"] {
            assert_eq!(score(text.as_bytes()), None, "{text:?}");
        }
    }

    #[test]
    fn an_mgyp_name_is_its_number() {
        assert_eq!(mgyp(b"MGYP000522683479"), Some(522_683_479));
        assert_eq!(mgyp(b"MGYP000000000000"), Some(0));
        assert_eq!(mgyp(b"MGYP999999999999"), Some(999_999_999_999));
        assert!(mgyp(b"MGYP999999999999").unwrap() < 1 << 40);

        for name in [
            &b"MGYP00052268347"[..],    // eleven digits
            b"MGYP0005226834799",       // thirteen
            b"MGYA000522683479",        // another prefix
            b"MGYP00052268347x",
            b"MGYP",
            b"7tm_1",
            b"",
        ] {
            assert_eq!(mgyp(name), None, "{:?}", String::from_utf8_lossy(name));
        }
    }
}
