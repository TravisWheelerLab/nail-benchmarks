//! Where recall's columns sit in a row.
//!
//! Everything about the file that is not those columns -- the preamble, the
//! blocks, the order, the trailer -- is [`super::frame`]'s, and the analyses
//! read a row through the frame. What is left here is the one thing that tells
//! this table from `runs.tbl`: a score column per tool rather than per run, so
//! the scores start after `pass` and there are as many of them as the file
//! declares tools.
//!
//! Nothing reads the score columns themselves, because no analysis here opens
//! one -- a summary asks only what the `pass` string says. One that wanted a
//! tool's score would add the accessor then, against this layout.

use super::Meta;
use super::frame::Layout;

/// Where recall's score columns start: after query, target and pass.
const SCORES: usize = 3;

/// Where this table's columns sit: query, target, pass, one per tool, inc,
/// dom.
pub fn layout(meta: &Meta) -> Layout {
    let tools = meta.tools().len();

    Layout {
        fields: SCORES + tools + 2,
        pass: 2,
        dom: SCORES + tools + 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::scores::FORMAT;
    use crate::scores::frame::Frame;

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

    /// A recall table, opened the way an analysis opens one.
    fn open(text: &str) -> anyhow::Result<Frame<&[u8]>> {
        let mut frame = Frame::new(text.as_bytes(), "test")?;

        anyhow::ensure!(frame.format() == FORMAT, "not a scores table");
        frame.layout(layout(&frame.meta));

        Ok(frame)
    }

    fn read(text: &str) -> anyhow::Result<Vec<String>> {
        let mut frame = open(text)?;
        let mut out = Vec::new();

        while frame.step()? {
            out.push(format!(
                "{} {} {} {} {}",
                String::from_utf8_lossy(frame.field(0)),
                String::from_utf8_lossy(frame.field(1)),
                String::from_utf8_lossy(frame.pass()),
                String::from_utf8_lossy(frame.field(SCORES)),
                frame.domain_count(),
            ));
        }

        Ok(out)
    }

    #[test]
    fn a_row_is_read_as_it_was_written() {
        assert_eq!(
            read(FILE).unwrap(),
            [
                "alpha MGYP000000000001 NH 25.0 1",
                "beta MGYP000000000002 Nh 11.0 0",
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
        let mut frame = open(FILE).unwrap();

        assert!(frame.step().unwrap());
        assert!(frame.passed(0));
        assert!(frame.passed(1));

        assert!(frame.step().unwrap());
        assert!(frame.passed(0));
        assert!(!frame.passed(1));
    }

    /// Every way the file can be wrong that would otherwise be read as an
    /// answer rather than as a fault.
    #[test]
    fn a_broken_file_is_refused() {
        let cases = [
            ("an older format", FILE.replace("#= format scores 2\n", "")),
            (
                "another table's format",
                FILE.replace("#= format scores 2", "#= format runs 1"),
            ),
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
            ("a row before any shard", FILE.replace("#= shard 1\n", "")),
        ];

        for (what, text) in cases {
            assert!(read(&text).is_err(), "{what} was read without complaint");
        }
    }
}
