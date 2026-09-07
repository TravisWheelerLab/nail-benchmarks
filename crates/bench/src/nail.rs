//! nail's hit table, for the one column `libsail`'s layout does not carry.
//!
//! ```text
//! |   0   |   1  |   2   |  3  |   4   |  5  |  6  |  7 |   8  |   9    |
//! # target  query  tstart tend  qstart qend  score bias evalue cellfrac
//! ```

use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::Context;

const N_FIELDS: usize = 10;

/// The fraction of the dynamic programming matrix nail computed for a hit.
pub struct CellFrac {
    pub query: String,
    pub target: String,
    pub cell_frac: f64,
}

/// Every row's cell fraction, in file order.
pub fn cell_fracs(path: impl AsRef<Path>) -> anyhow::Result<Vec<CellFrac>> {
    let path = path.as_ref();
    let reader = BufReader::new(
        std::fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?,
    );

    let mut out = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields = line.split_whitespace().collect::<Vec<_>>();
        if fields.len() < N_FIELDS {
            anyhow::bail!(
                "{} has a row of {} fields, nail writes {N_FIELDS}",
                path.display(),
                fields.len()
            );
        }

        out.push(CellFrac {
            target: fields[0].to_string(),
            query: fields[1].to_string(),
            cell_frac: fields[9].parse().with_context(|| {
                format!(
                    "unparseable cell fraction {:?} in {}",
                    fields[9],
                    path.display()
                )
            })?,
        });
    }

    Ok(out)
}
