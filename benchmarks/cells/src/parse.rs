//! Each form's nail table as the points the figure draws: one row per pair,
//! the cells of its full matrix against the fraction nail computed.
//!
//! What ran comes out of `ledger.tbl`, and which form a search covered out of
//! the set's `form` column, so a form is named by the recipe rather than by a
//! filename.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};
use clap::Parser;

use libsail::collection::Iterable;
use libsail::seq::fasta::IndexedFasta;
use libsail::seq::p7hmm::Hmm;
use libsail::tbl::Table;
use libsail::tbl::nail::NailRow;
use util::ledger::{self, Ledger};
use util::manifest;
use util::set::{Set, shape};

use crate::run::RUN_NAME;

#[derive(Parser)]
pub struct Args {
    /// Which label of paths.toml to read. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    /// Where the tables go. Defaults to what the label names
    #[arg(short, long, value_name = "dir")]
    out: Option<PathBuf>,
}

/// The table one form's points are written to.
pub fn table_name(form: &str) -> String {
    format!("cells.{form}.tbl")
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    let set = Set::load_as(&paths.set, &shape::FORMS)?;

    let ran = Ledger::load(&paths.run)?;
    ledger::warn(ran.failed(), "form(s)");

    let out = args.out.unwrap_or_else(|| paths.analysis.clone());
    std::fs::create_dir_all(&out).with_context(|| format!("failed to create {}", out.display()))?;

    let results = paths.run.join("results");
    let mut query_lens: HashMap<PathBuf, HashMap<String, usize>> = HashMap::new();
    let mut parsed = 0usize;

    for row in ran.runs().filter(|row| row.name == RUN_NAME) {
        ensure!(
            row.tool == "nail",
            "run {RUN_NAME:?} was produced by {}, and only nail reports cell fractions",
            row.tool
        );

        let unit = set
            .units()
            .find(|u| u.name() == row.shard)
            .with_context(|| format!("the set has no unit {:?}", row.shard))?;
        let form = unit.need("form")?;

        // every form is searched by the same profiles, so they are read once
        let hmm = unit.query_hmm()?;
        if !query_lens.contains_key(&hmm) {
            query_lens.insert(hmm.clone(), model_lengths(&hmm)?);
        }
        let query_lens = &query_lens[&hmm];
        let target_lens = seq_lengths(&unit.target()?)?;

        let table = manifest::table_path(&results, RUN_NAME, &row.shard);
        let hits = Table::<NailRow>::open(&table)
            .with_context(|| format!("failed to read {}", table.display()))?;

        let mut points = toil::Table::new(toil::Schema::new(["cells", "fraction"]));
        for h in hits.iter() {
            let qlen = query_lens
                .get(&h.query)
                .with_context(|| format!("no query len for: {}", h.query))?;
            let tlen = target_lens
                .get(&h.target)
                .with_context(|| format!("no target len for: {}", h.target))?;

            points.row([
                toil::Cell::from(qlen * tlen),
                toil::Cell::from(h.cell_fraction),
            ]);
        }

        let path = out.join(table_name(form));
        points
            .write(&path)
            .with_context(|| format!("failed to write {}", path.display()))?;
        println!("wrote {}", path.display());
        parsed += 1;
    }

    if parsed == 0 {
        bail!("no finished {RUN_NAME:?} runs in {}", paths.run.display());
    }

    Ok(())
}

/// Each profile's model length, by name.
//
// read out of the file rather than shelled out to hmmstat, which is not a
// tool this benchmark builds
fn model_lengths(path: &Path) -> anyhow::Result<HashMap<String, usize>> {
    Ok(Hmm::open(path)
        .with_context(|| format!("failed to parse {}", path.display()))?
        .iter()
        .map(|model| {
            (
                String::from_utf8_lossy(&model.header.name).into_owned(),
                model.header.leng,
            )
        })
        .collect())
}

/// Each sequence's length, by name.
fn seq_lengths(path: &Path) -> anyhow::Result<HashMap<String, usize>> {
    let fa =
        IndexedFasta::open(path).with_context(|| format!("failed to open {}", path.display()))?;

    fa.iter()
        .map(|rec| Ok((rec.name_str()?.to_string(), rec.seq.len())))
        .collect()
}
