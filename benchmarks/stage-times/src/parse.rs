//! Reads each unit's stage tree into `stages.tbl`, and its counts into
//! `counts.tbl`.
//!
//! nail prints the two branches of the tree as wall clock, and each leaf as the
//! sum of that stage's time over every pair, across every thread. So a leaf's
//! percentage is its share of the thread-seconds its branch spent, and that
//! share of the branch's wall is the wall the leaf gets here. It is exact when
//! every thread is busy for the whole branch; idle time at the tail of a
//! branch is spread over the leaves in proportion rather than shown.

use std::path::PathBuf;

use anyhow::{Context, bail, ensure};
use clap::Parser;

use util::ledger::{self, Ledger};
use util::search::Dirs;

use crate::run::{RUN_NAME, stats_path};

#[derive(Parser)]
pub struct Args {
    /// Which label of paths.toml to read. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,
}

/// One line of the runtime tree.
struct Timed {
    /// The stage this line is nested under; empty for the total.
    parent: String,
    stage: String,
    secs: f64,
    /// Its share of what its parent printed, in [0, 1].
    share: f64,
}

/// One line of the counts block.
struct Counted {
    quantity: String,
    value: u64,
}

pub fn main(_: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    let dirs = Dirs::new(&paths.run, &paths.tmp);
    let units = units(&dirs)?;

    std::fs::create_dir_all(&paths.analysis)
        .with_context(|| format!("failed to create {}", paths.analysis.display()))?;

    let mut stages = toil::Table::new(toil::Schema::new([
        toil::Column::new("unit"),
        toil::Column::new("parent"),
        toil::Column::new("stage"),
        toil::Column::new("secs").fixed(2),
        toil::Column::new("share").fixed(4),
        toil::Column::new("wall").fixed(2),
    ]));
    let mut counts = toil::Table::new(toil::Schema::new([
        toil::Column::new("unit"),
        toil::Column::new("quantity"),
        toil::Column::new("value"),
    ]));

    for (unit, wall_s, path) in units {
        let text = std::fs::read_to_string(&path)
            .with_context(|| format!("failed to read {}", path.display()))?;
        let (timed, counted) =
            read(&text).with_context(|| format!("failed to parse {}", path.display()))?;

        // what michi timed for the whole command, beside what nail
        // printed for itself
        stages.meta("michi", [unit.clone(), format!("{wall_s:.2}")]);

        // the total and the branches are wall as nail printed them;
        // a leaf is its share of its branch's
        let mut wall_of = std::collections::HashMap::new();
        for t in &timed {
            let wall = match wall_of.get(&t.parent) {
                Some(branch) if t.parent != "total" => t.share * branch,
                _ => t.secs,
            };
            wall_of.insert(t.stage.clone(), wall);
            stages.row([
                toil::Cell::from(&unit),
                t.parent.as_str().into(),
                t.stage.as_str().into(),
                t.secs.into(),
                t.share.into(),
                wall.into(),
            ]);
        }

        for c in counted {
            counts.row([
                toil::Cell::from(&unit),
                c.quantity.as_str().into(),
                c.value.into(),
            ]);
        }
    }

    for (table, name) in [(&stages, "stages.tbl"), (&counts, "counts.tbl")] {
        let path = paths.analysis.join(name);
        table
            .write(&path)
            .with_context(|| format!("failed to write {}", path.display()))?;
        println!("wrote {}", path.display());
    }
    Ok(())
}

/// The units the run finished, with michi's wall clock for each and where
/// nail's stdout went, in the order the pipeline declared them.
fn units(dirs: &Dirs) -> anyhow::Result<Vec<(String, f64, PathBuf)>> {
    let ran = Ledger::load(&dirs.root)?;
    ledger::warn(ran.failed(), "unit(s)");

    let mut units = Vec::new();
    for row in ran.runs().filter(|row| row.name == RUN_NAME) {
        ensure!(
            row.tool == "nail",
            "run {RUN_NAME:?} was produced by {}",
            row.tool
        );
        ensure!(
            !row.shard.is_empty(),
            "run {RUN_NAME:?} has no shard saying which unit it was"
        );

        let wall_s = row.wall_s.context("no wall clock")?;
        units.push((row.shard.clone(), wall_s, stats_path(dirs, &row.shard)));
    }

    if units.is_empty() {
        bail!("no finished {RUN_NAME:?} runs in {}", dirs.root.display());
    }
    Ok(units)
}

/// The runtime tree and the counts block out of what nail printed.
//
// everything above "summary statistics:" is progress output, some of
// it overwritten in place with cursor escapes, and is skipped. the
// tree is read by indent: the total sits flush left, a branch is
// indented under it, a leaf under a branch
fn read(text: &str) -> anyhow::Result<(Vec<Timed>, Vec<Counted>)> {
    let body = text
        .split_once("summary statistics:")
        .map(|(_, rest)| rest)
        .context("no summary statistics block")?;
    let (counts, runtime) = body.split_once("runtime:").context("no runtime block")?;

    let counted = counts
        .lines()
        .filter_map(|line| line.split_once(':'))
        .map(|(label, value)| {
            let value = value.trim().replace(',', "");
            Ok(Counted {
                quantity: snake(label),
                value: value
                    .parse()
                    .with_context(|| format!("bad count {value:?}"))?,
            })
        })
        .collect::<anyhow::Result<Vec<_>>>()?;

    let mut timed = Vec::new();
    let mut branch = String::new();
    for (i, line) in format!("runtime:{runtime}").lines().enumerate() {
        let Some((label, value)) = line.rsplit_once(':') else {
            continue;
        };
        let indent = line.len() - line.trim_start().len();
        let stage = snake(label);
        let parent = match (i, indent) {
            (0, _) => String::new(),
            (_, i) if i < 4 => {
                branch = stage.clone();
                "total".to_string()
            }
            _ => branch.clone(),
        };

        let (secs, pct) = value.trim().split_once('s').context("no seconds")?;
        let pct = pct
            .trim()
            .trim_matches(|c| c == '(' || c == ')' || c == '%')
            .trim();

        timed.push(Timed {
            parent,
            stage,
            secs: secs
                .parse()
                .with_context(|| format!("bad seconds {secs:?}"))?,
            share: pct
                .parse::<f64>()
                .with_context(|| format!("bad percent {pct:?}"))?
                / 100.0,
        });
    }

    ensure!(!timed.is_empty(), "empty runtime block");
    Ok((timed, counted))
}

/// A tree label as a cell: the tree glyphs and punctuation gone, words joined
/// by underscores, and the parenthetical dropped where it only says which tool.
fn snake(label: &str) -> String {
    let label = label
        .trim_start_matches(|c: char| c.is_whitespace() || "├└─".contains(c))
        .trim();
    let label = match label {
        "runtime" => "total",
        "seeding (mmseqs)" => "seeding",
        "[misc.]" => "misc",
        other => other,
    };
    label
        .to_lowercase()
        .split(|c: char| !c.is_alphanumeric())
        .filter(|w| !w.is_empty())
        .collect::<Vec<_>>()
        .join("_")
}

#[cfg(test)]
mod tests {
    use super::*;

    const OUT: &str = "\
indexing query database... done (0.01s)
summary statistics:
 ├─ queries:                                   200
 ├─ targets:                                20,000
 └─ backward DP cells computed:         38,674,118

runtime: 4.00s (100.00%)
 └─ seeding (mmseqs):   1.89s (47.13%)
     ├─ prefilter:     1.21s (63.95%)
     └─ [misc.]:       0.23s (12.33%)
 └─ alignment:          1.96s (49.02%)
     ├─ cloud search:       7.71s (55.32%)
     ├─ output (mutex):     0.00s ( 0.00%)
     └─ [misc.]:            0.23s ( 1.62%)
";

    #[test]
    fn reads_the_tree_by_indent() {
        let (timed, counted) = read(OUT).unwrap();

        // shares back to the two places nail printed, so the
        // comparison is on text rather than on floats
        let rows: Vec<(&str, &str, f64, String)> = timed
            .iter()
            .map(|t| {
                (
                    t.parent.as_str(),
                    t.stage.as_str(),
                    t.secs,
                    format!("{:.2}", t.share * 100.0),
                )
            })
            .collect();
        assert_eq!(
            rows,
            [
                ("", "total", 4.00, "100.00".into()),
                ("total", "seeding", 1.89, "47.13".into()),
                ("seeding", "prefilter", 1.21, "63.95".into()),
                ("seeding", "misc", 0.23, "12.33".into()),
                ("total", "alignment", 1.96, "49.02".into()),
                ("alignment", "cloud_search", 7.71, "55.32".into()),
                ("alignment", "output_mutex", 0.00, "0.00".into()),
                ("alignment", "misc", 0.23, "1.62".into()),
            ]
        );

        let counts: Vec<(&str, u64)> = counted
            .iter()
            .map(|c| (c.quantity.as_str(), c.value))
            .collect();
        assert_eq!(
            counts,
            [
                ("queries", 200),
                ("targets", 20000),
                ("backward_dp_cells_computed", 38674118),
            ]
        );
    }
}
