//! The padded table every analysis in this repo writes.
//!
//! A `#` header naming the columns, a `#` rule of dashes under it, then one row
//! per line with the cells padded so they sit under their names. Above all of
//! that, an optional block of `#` metadata: what was searched, what a fraction
//! is a fraction of, whatever the figures want as a reference line.
//!
//! It is a format for reading, not for round-tripping. A reader that wants its
//! own table back writes the parser for it -- `#=` metadata lines and a known
//! column order are what make that possible, and neither is this module's
//! business.

use std::path::Path;

use anyhow::Context;

pub struct Table<'a> {
    /// Lines above the header, each already `#`-prefixed and newline-ended.
    /// Empty for a table that has nothing to say about itself.
    pub meta: &'a str,
    pub headers: &'a [String],
    /// One per row, each as wide as `headers`.
    pub rows: &'a [Vec<String>],
    /// Leave the last column unpadded.
    ///
    /// For a column whose width is a property of the row rather than of the
    /// table -- a comma-joined list of domain scores, an argv. Padding it to
    /// the widest row would put most of a line's bytes into trailing spaces.
    pub ragged_last: bool,
}

pub fn write(path: &Path, table: Table<'_>) -> anyhow::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }

    std::fs::write(path, render(table))
        .with_context(|| format!("failed to write {}", path.display()))
}

/// The table as text, for a caller that has somewhere else to put it.
pub fn render(table: Table<'_>) -> String {
    let Table {
        meta,
        headers,
        rows,
        ragged_last,
    } = table;

    let widths: Vec<usize> = headers
        .iter()
        .enumerate()
        .map(|(i, h)| {
            rows.iter()
                .filter_map(|r| r.get(i))
                .map(String::len)
                .chain(std::iter::once(h.len()))
                .max()
                .unwrap_or(0)
        })
        .collect();

    // no columns means no last column to leave ragged
    let last = ragged_last.then(|| headers.len().checked_sub(1)).flatten();

    let mut out = meta.to_string();

    out.push('#');
    for (i, (h, &w)) in headers.iter().zip(&widths).enumerate() {
        match last == Some(i) {
            true => out.push_str(&format!(" {h}")),
            false => out.push_str(&format!(" {h:<w$}")),
        }
    }

    out.push_str("\n#");
    for (i, &w) in widths.iter().enumerate() {
        // the rule under a ragged column is as wide as its name, since there is
        // no column width for it to be as wide as
        let w = match last == Some(i) {
            true => headers[i].len(),
            false => w,
        };
        out.push_str(&format!(" {}", "-".repeat(w)));
    }
    out.push('\n');

    for row in rows {
        // the two characters the `# ` takes on a header line, so the cells sit
        // under their names rather than beside them
        out.push_str("  ");
        for (i, (c, &w)) in row.iter().zip(&widths).enumerate() {
            if i > 0 {
                out.push(' ');
            }
            match last == Some(i) {
                true => out.push_str(c),
                false => out.push_str(&format!("{c:<w$}")),
            }
        }
        out.push('\n');
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(xs: &[&str]) -> Vec<String> {
        xs.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn cells_sit_under_their_names() {
        let headers = strings(&["name", "n"]);
        let rows = vec![strings(&["a-long-one", "1"]), strings(&["b", "22"])];

        let text = render(Table {
            meta: "",
            headers: &headers,
            rows: &rows,
            ragged_last: false,
        });

        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "# name       n ");
        assert_eq!(lines[1], "# ---------- --");
        assert_eq!(lines[2], "  a-long-one 1 ");
        assert_eq!(lines[3], "  b          22");

        // the first column starts at the same offset on every line, which is
        // the whole point of the two leading spaces on a row
        for line in &lines {
            assert!(
                matches!(&line[..2], "# " | "  "),
                "line does not start its first column at offset 2: {line:?}"
            );
        }
    }

    #[test]
    fn a_ragged_last_column_is_not_padded() {
        let headers = strings(&["name", "doms"]);
        let rows = vec![strings(&["a", "1.0,2.0,3.0"]), strings(&["b", "4.0"])];

        let text = render(Table {
            meta: "",
            headers: &headers,
            rows: &rows,
            ragged_last: true,
        });

        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "# name doms");
        assert_eq!(lines[1], "# ---- ----");
        assert_eq!(lines[3], "  b    4.0");
    }

    #[test]
    fn meta_goes_above_the_header() {
        let headers = strings(&["a"]);
        let text = render(Table {
            meta: "# query 200\n#\n",
            headers: &headers,
            rows: &[],
            ragged_last: false,
        });

        assert!(text.starts_with("# query 200\n#\n# a\n# -\n"));
    }
}
