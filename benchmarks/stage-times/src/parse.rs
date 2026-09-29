//! Reads each unit's stage tree into `stages.tbl`, and its counts into
//! `counts.tbl`.
//!
//! nail prints every line of the tree as wall clock. Under `alignment`, which
//! runs on every thread, it also prints each stage's cpu seconds, the sum of
//! that stage's time over every pair across every thread; the wall there is
//! nail's own share-of-cpu times branch wall. Under `align` in prog mode it
//! prints one line per round.

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
#[derive(Debug, PartialEq)]
struct Timed {
    /// The line this one is nested under; empty for the total.
    parent: String,
    stage: String,
    wall: f64,
    /// Its share of its parent's wall, in [0, 1].
    share: f64,
    /// Seconds on cpu, where nail prints them.
    cpu: Option<f64>,
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
        toil::Column::new("wall").fixed(2),
        toil::Column::new("share").fixed(4),
        toil::Column::new("cpu").fixed(2),
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

        for t in &timed {
            stages.row([
                toil::Cell::from(&unit),
                t.parent.as_str().into(),
                t.stage.as_str().into(),
                t.wall.into(),
                t.share.into(),
                t.cpu.map(toil::Cell::from).unwrap_or(toil::Cell::from("")),
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
// it overwritten in place with cursor escapes, and is skipped
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

    let timed = tree(&format!("runtime:{runtime}"))?;
    ensure!(!timed.is_empty(), "empty runtime block");
    Ok((timed, counted))
}

/// The runtime tree, read by indent.
//
// the total sits flush left, a branch is indented under it, a leaf
// under a branch, and a round under the align leaf. what the line is
// nested under is the last line seen one level up. the alignment
// block's column header carries no number and is skipped with the
// blank lines
fn tree(runtime: &str) -> anyhow::Result<Vec<Timed>> {
    // the stage name at each depth, so a line's parent is the name
    // one level up
    let mut above: Vec<String> = Vec::new();
    // the wall at each depth, for a line's share of its parent
    let mut walls: Vec<f64> = Vec::new();
    let mut timed = Vec::new();

    for line in runtime.lines() {
        let Some(row) = parse_line(line)? else {
            continue;
        };

        let depth = row.depth;
        above.truncate(depth);
        walls.truncate(depth);

        let parent = above.last().cloned().unwrap_or_default();
        let share = match (row.share, walls.last()) {
            (Some(share), _) => share,
            (None, Some(parent_wall)) if *parent_wall > 0.0 => row.wall / parent_wall,
            (None, _) => 1.0,
        };

        above.push(row.stage.clone());
        walls.push(row.wall);

        timed.push(Timed {
            parent,
            stage: row.stage,
            wall: row.wall,
            share,
            cpu: row.cpu,
        });
    }

    Ok(timed)
}

/// One line of the tree with its numbers read, before it knows its parent.
struct Line {
    depth: usize,
    stage: String,
    wall: f64,
    share: Option<f64>,
    cpu: Option<f64>,
}

/// Read one line, or `None` for a line that is not a stage.
//
// the label ends at the first word that is a number of seconds. after
// it come an optional percentage, in parentheses or bare, an optional
// second number of seconds which is cpu, and an optional [note]
fn parse_line(line: &str) -> anyhow::Result<Option<Line>> {
    let glyphs: &[char] = &[' ', '│', '├', '└', '─'];
    let text = line.trim_start_matches(glyphs);
    if text.is_empty() {
        return Ok(None);
    }

    // four columns of indent per level of the tree; the glyphs are
    // several bytes each, so this counts chars
    let depth = (line.chars().count() - text.chars().count()) / 4;

    let words: Vec<&str> = text.split_whitespace().collect();
    let Some(at) = words.iter().position(|w| seconds(w).is_some()) else {
        return Ok(None);
    };
    if at == 0 {
        return Ok(None);
    }

    let label = words[..at].join(" ");
    let wall = seconds(words[at]).unwrap();

    // the note after a branch: "[6 iterations]", or "[4 threads, cpu
    // 147.31s, 99.5% busy]", where the cpu is worth keeping. it comes
    // after the seconds, so "[misc.]" as a label is not one
    let rest = &words[at + 1..];
    let note_at = rest
        .iter()
        .position(|w| w.starts_with('['))
        .unwrap_or(rest.len());
    let (rest, note) = rest.split_at(note_at);

    let mut share = None;
    let mut cpu = note
        .iter()
        .map(|w| w.trim_matches(|c| c == '[' || c == ']' || c == ','))
        // "cpu 7,842.90s," keeps its inner comma and loses the last
        .skip_while(|w| *w != "cpu")
        .nth(1)
        .and_then(seconds);
    for w in rest {
        let w = w.trim_matches(|c| c == '(' || c == ')');
        if let Some(pct) = w.strip_suffix('%') {
            share = Some(
                pct.parse::<f64>()
                    .with_context(|| format!("bad percent {w:?}"))?
                    / 100.0,
            );
        } else if let Some(s) = seconds(w) {
            cpu = Some(s);
        }
    }

    Ok(Some(Line {
        depth,
        stage: snake(label.trim_end_matches(':')),
        wall,
        share,
        cpu,
    }))
}

/// A word like `12.34s` or `7,842.90s` as seconds.
fn seconds(word: &str) -> Option<f64> {
    word.strip_suffix('s')?.replace(',', "").parse().ok()
}

/// A tree label as a cell: punctuation gone, words joined by underscores.
fn snake(label: &str) -> String {
    let label = match label.trim() {
        "runtime" => "total",
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

runtime: 566.03s (100.00%)
 └─ setup:       3.91s ( 0.69%)
     ├─ query index:      0.79s (20.13%)
     └─ [misc.]:          0.00s ( 0.02%)
 └─ seeding:   525.11s (92.77%)
     ├─ prefilter:   392.68s (74.78%)
     ├─ align:       128.71s (24.51%)   [2 iterations]
     │   ├─ 1:        64.85s
     │   └─ 2:        63.86s
     └─ [misc.]:       0.02s ( 0.00%)
 └─ alignment:  37.01s ( 6.54%)   [4 threads, cpu 1,147.31s, 99.5% busy]
     │                       wall        %      cpu
     │                     ------   ------   ------
     ├─ memory init         0.18s    0.48%    0.71s
     ├─ cloud search       16.65s   44.98%   1,066.26s
     └─ [misc.]             0.49s    1.32%    1.95s
 └─ [misc.]:     0.00s ( 0.00%)
";

    fn row(parent: &str, stage: &str, wall: f64, share: f64, cpu: Option<f64>) -> Timed {
        Timed {
            parent: parent.into(),
            stage: stage.into(),
            wall,
            share,
            cpu,
        }
    }

    #[test]
    fn reads_the_tree_by_indent() {
        let (timed, counted) = read(OUT).unwrap();

        // shares back to the two places nail printed, and the rounds'
        // shares to four, so the comparison is on text rather than on
        // floats
        let rounded: Vec<Timed> = timed
            .into_iter()
            .map(|t| Timed {
                share: (t.share * 10000.0).round() / 10000.0,
                ..t
            })
            .collect();

        assert_eq!(
            rounded,
            [
                row("", "total", 566.03, 1.0, None),
                row("total", "setup", 3.91, 0.0069, None),
                row("setup", "query_index", 0.79, 0.2013, None),
                row("setup", "misc", 0.00, 0.0002, None),
                row("total", "seeding", 525.11, 0.9277, None),
                row("seeding", "prefilter", 392.68, 0.7478, None),
                row("seeding", "align", 128.71, 0.2451, None),
                row("align", "1", 64.85, 0.5038, None),
                row("align", "2", 63.86, 0.4962, None),
                row("seeding", "misc", 0.02, 0.0, None),
                row("total", "alignment", 37.01, 0.0654, Some(1147.31)),
                row("alignment", "memory_init", 0.18, 0.0048, Some(0.71)),
                row("alignment", "cloud_search", 16.65, 0.4498, Some(1066.26)),
                row("alignment", "misc", 0.49, 0.0132, Some(1.95)),
                row("total", "misc", 0.00, 0.0, None),
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
