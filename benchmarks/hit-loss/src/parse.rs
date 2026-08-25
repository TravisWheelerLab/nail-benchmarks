//! Turns a finished run into a funnel: how many of hmmer's hits survive each
//! checkpoint nail's pipeline offers a per-pair view of.
//!
//! Only two checkpoints are directly observable without instrumenting nail
//! itself: whether a pair got a seed, and whether it ended up in nail's
//! .tbl. Everything cloud search and alignment do to a seeded pair between
//! those two points is invisible here and collapses into one bucket, "seeded
//! but unreported" -- see run.rs for why nail's -E has to be cranked up for
//! that bucket to mean what it says.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, bail};
use clap::Parser;

use bioio::tbl::{HitTable, HmmerTable, NailTable};

type Pair = (String, String);

#[derive(Parser, Debug)]
pub struct Args {
    /// A benchmark directory `hit-loss run` has finished on.
    bench_dir: PathBuf,

    /// The per-family cutoffs mgnify learned. Defaults to its committed ones.
    #[arg(long, value_name = "cutoffs.txt")]
    cutoffs: Option<PathBuf>,

    /// Which decoy to cut at. The cutoffs file holds each family's five
    /// best-scoring decoys, so this admits at most `c` false positives per
    /// family.
    #[arg(short = 'c', default_value_t = 2, value_name = "N")]
    c: usize,

    #[arg(short, long, value_name = "funnel.tbl")]
    out: Option<PathBuf>,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let cutoffs_path = args
        .cutoffs
        .unwrap_or_else(|| crate::dir().join("../mgnify/cutoffs.txt"));
    let out = args
        .out
        .unwrap_or_else(|| args.bench_dir.join("funnel.tbl"));

    let cutoffs = cutoffs(&cutoffs_path, args.c)
        .with_context(|| format!("failed to read {}", cutoffs_path.display()))?;

    let hmmer = HitTable::from_path::<_, HmmerTable>(args.bench_dir.join("truth/hmmer.tbl"))
        .context("failed to read truth/hmmer.tbl; has `hit-loss run` finished?")?;

    let seeded = pairs(args.bench_dir.join("seeds/seeds"))
        .context("failed to read seeds/seeds; has `hit-loss run` finished?")?;

    let reported: HashSet<Pair> =
        HitTable::from_path::<_, NailTable>(args.bench_dir.join("nail.tbl"))
            .context("failed to read nail.tbl; has `hit-loss run` finished?")?
            .hits
            .into_iter()
            .map(|h| (h.query, h.target))
            .collect();

    let mut truth = 0usize;
    let mut lost_seed = 0usize;
    let mut lost_cloud_align = 0usize;
    let mut reached = 0usize;

    for hit in &hmmer.hits {
        let Some(&cutoff) = cutoffs.get(&hit.query) else {
            continue;
        };
        if hit.score < cutoff {
            continue;
        }

        truth += 1;

        let pair = (hit.query.clone(), hit.target.clone());
        if !seeded.contains(&pair) {
            lost_seed += 1;
        } else if !reported.contains(&pair) {
            lost_cloud_align += 1;
        } else {
            reached += 1;
        }
    }

    ensure_nonzero(truth)?;

    write(&out, truth, lost_seed, lost_cloud_align, reached)?;

    println!(
        "wrote {} ({truth} true hits: {lost_seed} lost at seeding, \
         {lost_cloud_align} lost in cloud/align, {reached} reported)",
        out.display()
    );

    Ok(())
}

fn ensure_nonzero(truth: usize) -> anyhow::Result<()> {
    if truth == 0 {
        bail!("no hmmer hit cleared its family's cutoff; there is nothing to measure against");
    }
    Ok(())
}

/// The (query, target) pairs `--seeds-out` wrote: nail's own prf/seq column
/// order, whitespace-separated, no header.
fn pairs(path: PathBuf) -> anyhow::Result<HashSet<Pair>> {
    let file = File::open(&path).with_context(|| format!("failed to open {}", path.display()))?;

    BufReader::new(file)
        .lines()
        .map(|line| {
            let line = line?;
            let mut fields = line.split_whitespace();
            let query = fields.next().context("a seed row has no query")?;
            let target = fields.next().context("a seed row has no target")?;
            Ok((query.to_string(), target.to_string()))
        })
        .collect()
}

type Cutoffs = std::collections::HashMap<String, f32>;

/// One score per family, out of the decoys mgnify scored it against.
///
/// A family is kept only when both nail and mmseqs got a nonzero cutoff,
/// which is mgnify's rule and is copied rather than improved so the numbers
/// stay comparable to its own.
fn cutoffs(path: &Path, c: usize) -> anyhow::Result<Cutoffs> {
    let reader = BufReader::new(File::open(path)?);
    let mut out = Cutoffs::new();

    for line in reader.lines() {
        let line = line?;

        let Some((family, rest)) = line.split_once(',') else {
            continue;
        };

        let groups: Vec<(&str, Vec<f32>)> = rest
            .split("),(")
            .map(|g| g.trim_matches(|c| c == '(' || c == ')'))
            .map(|g| {
                let mut it = g.split(',');
                let tool = it.next().unwrap_or_default();
                let mut nums: Vec<f32> = it.filter_map(|x| x.parse().ok()).collect();
                // the last number is how many decoys there were, not a score
                nums.pop();
                (tool, nums)
            })
            .collect();

        let nail = groups.iter().find(|(t, _)| *t == "nail");
        let mmseqs = groups.iter().find(|(t, _)| *t == "mmseqs");

        let (Some((_, nail)), Some((_, mmseqs))) = (nail, mmseqs) else {
            continue;
        };

        if let (Some(&n), Some(&m)) = (nail.get(c), mmseqs.get(c))
            && n > 0.0
            && m > 0.0
        {
            out.insert(family.to_string(), n);
        }
    }

    if out.is_empty() {
        bail!("no usable cutoffs at index {c}");
    }

    Ok(out)
}

fn write(
    path: &Path,
    truth: usize,
    lost_seed: usize,
    lost_cloud_align: usize,
    reached: usize,
) -> anyhow::Result<()> {
    let rows = [
        ("truth", truth, truth),
        ("lost_seed", lost_seed, truth - lost_seed),
        (
            "lost_cloud_align",
            lost_cloud_align,
            truth - lost_seed - lost_cloud_align,
        ),
        ("reported", reached, reached),
    ];

    let headers = ["stage", "n", "sens"];

    let cells: Vec<[String; 3]> = rows
        .iter()
        .map(|&(stage, n, remaining)| {
            [
                stage.to_string(),
                n.to_string(),
                format!("{:.4}", remaining as f64 / truth as f64),
            ]
        })
        .collect();

    let widths: Vec<usize> = headers
        .iter()
        .enumerate()
        .map(|(i, h)| {
            cells
                .iter()
                .map(|c| c[i].len())
                .chain(std::iter::once(h.len()))
                .max()
                .unwrap_or(0)
        })
        .collect();

    let mut out = String::new();
    let pad = |text: &str, width: usize| format!("{text:<width$}");

    out.push_str("# ");
    for (h, &w) in headers.iter().zip(&widths) {
        out.push_str(&pad(h, w));
        out.push(' ');
    }
    out.push_str("\n# ");
    for &w in &widths {
        out.push_str(&"-".repeat(w));
        out.push(' ');
    }
    out.push('\n');

    for row in &cells {
        out.push_str("  ");
        for (c, &w) in row.iter().zip(&widths) {
            out.push_str(&pad(c, w));
            out.push(' ');
        }
        out.push('\n');
    }

    std::fs::write(path, out).with_context(|| format!("failed to write {}", path.display()))
}
