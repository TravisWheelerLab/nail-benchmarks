//! Cuts a query source and a target source into a set the benchmarks can
//! search.
//!
//! A build produces one directory under `store/`, holding the files and a
//! `set.tbl` describing them. Nothing downstream reads the directory layout:
//! a pipeline loads the manifest and gets the query, the target and whatever
//! the recipe wrote down about each unit, so a search does not learn which of
//! the two recipes ran.
//!
//! The two recipes differ in how they cut, not in what they cut. `fixed` deals
//! one query set against target shards of equal size, which is what a question
//! about recall needs: the shards are units of work rather than a variable.
//! `ladder` nests rungs on both axes, each a prefix of the one above, which is
//! what a question about scaling needs: every rung is a measurement.
//!
//! Neither knows what it is cutting, nor where the result goes. A label in
//! `paths.toml` names the sources, the size and the destination, and this
//! stamps the shape into the `set.tbl` it leaves there. So the benchmarks read
//! a contract rather than a location, and cutting something other than Pfam is
//! editing a line.
//!
//! What they share is here: subsetting the query, and building the mmseqs
//! profile db from what came out.
//!
//! The stockholm alignments and the mmseqs profile db are always built, even
//! for the pipelines that don't search mmseqs. They are what an mmseqs column
//! costs, and a query set that can't answer for one of the tools isn't one set.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use anyhow::{Context, bail, ensure};
use clap::Parser;
use serde::Deserialize;
use libsail::collection::{Aggregate, Indexable, Iterable};
use libsail::format::{self, Format};
use libsail::index::{self, Index};
use libsail::seq::fasta::{DEFAULT_LINE_WIDTH, IndexedFasta};
use libsail::seq::p7hmm::IndexedHmm;
use rand::SeedableRng;
use rand::rngs::StdRng;

use michi::{Closure, Cmd as PCmd, PipelineBuilder, Progress, Step};
mod profmark;

use util::cut;
use util::paths;
use util::set::{self, Set};
use util::tools::mmseqs;


/// Where a set keeps things, relative to its own root.
//
// these spell both the paths a build writes and the cells the manifest
// carries, so the two cannot drift apart
const QUERIES: &str = "queries";
const TARGETS: &str = "targets";
const QUERY_HMM: &str = "query.hmm";
const QUERY_STO: &str = "query.sto";
const QUERY_DB: &str = "queryDB/queryDB";

// --------------------------------------------------------------------- cli

/// One dataset this can make: where it comes from, how much of it, where it
/// goes.
///
/// Tagged by `shape`, so a label says which recipe it is and there is one place
/// that is named. `deny_unknown_fields` is what makes a `fixed` entry carrying
/// `target-rungs` a parse error naming the line, rather than a key silently
/// ignored an hour into a draw.
#[derive(Deserialize, Debug)]
#[serde(
    tag = "shape",
    rename_all = "kebab-case",
    rename_all_fields = "kebab-case",
    deny_unknown_fields
)]
enum Recipe {
    /// One query set against target shards of equal size.
    Fixed {
        queries: PathBuf,
        alignments: PathBuf,
        targets: PathBuf,
        out: PathBuf,
        shards: usize,
        /// A cap on the sequences drawn. All of them when absent.
        #[serde(default)]
        seqs: Option<usize>,
        /// A cap on the families cut. All of them when absent.
        #[serde(default)]
        fams: Option<usize>,
        /// Write every sequence backwards.
        ///
        /// The same deal, not a pass over a finished set: the draw, the
        /// sharding and the seed are untouched, so a reversed set holds the
        /// reverse of exactly what its forward twin holds. Reversing keeps a
        /// sequence's composition and destroys its homology, which is what a
        /// calibration searches against.
        #[serde(default)]
        reversed: bool,
        #[serde(default = "default_seed")]
        seed: u64,
    },
    /// One query against one target, paired rather than crossed.
    ///
    /// The two sources are separate directories of fasta, so a sequence is
    /// never on both sides: pair `i` is the `i`th file of each, in name order.
    /// Nothing is drawn -- these are pairs somebody chose, and the recipe
    /// selects and places them.
    Pairs {
        queries: PathBuf,
        targets: PathBuf,
        out: PathBuf,
        /// How many of the pairs to take, in name order. All of them when
        /// absent.
        #[serde(default)]
        pairs: Option<usize>,
    },
    /// One search over a profmark split, with a truth table beside it.
    ///
    /// The split is expensive and depends only on the alignments and the split
    /// parameters, so it is drawn once into `split` and reused. The set itself
    /// is one unit: every query against one target file, with `truth.tbl`
    /// saying which pair is which and at what identity.
    Profmark {
        /// The alignments the split is drawn from, Pfam's SEED.
        alignments: PathBuf,
        /// The sequences the true targets are hidden among, Swissprot.
        decoys: PathBuf,
        /// Where the train/test split lives, built if it is not there.
        split: PathBuf,
        out: PathBuf,
        /// A cap on the benchmark pairs. All that survive filtering when
        /// absent; decoys go in on top of this.
        #[serde(default)]
        pairs: Option<usize>,
        #[serde(default = "default_seed")]
        seed: u64,
        /// Maximum identity between the train and test halves.
        #[serde(default = "default_train_test_id")]
        train_test_id: f64,
        #[serde(default = "default_min_test")]
        min_test: usize,
        #[serde(default = "default_max_test")]
        max_test: usize,
        /// Threads for hmmbuild.
        #[serde(default = "default_threads")]
        threads: usize,
    },
    /// Nested rungs on both axes, each a prefix of the one above.
    Ladder {
        queries: PathBuf,
        alignments: PathBuf,
        targets: PathBuf,
        out: PathBuf,
        query_rungs: Vec<usize>,
        target_rungs: Vec<usize>,
        #[serde(default = "default_seed")]
        seed: u64,
    },
}

fn default_train_test_id() -> f64 {
    0.30
}

fn default_min_test() -> usize {
    10
}

fn default_max_test() -> usize {
    30
}

fn default_threads() -> usize {
    8
}

fn default_seed() -> u64 {
    67779
}

impl Recipe {
    /// Where this recipe puts what it makes.
    fn out(&self) -> &Path {
        match self {
            Recipe::Fixed { out, .. }
            | Recipe::Pairs { out, .. }
            | Recipe::Ladder { out, .. }
            | Recipe::Profmark { out, .. } => out,
        }
    }

    fn shape(&self) -> &'static str {
        match self {
            Recipe::Fixed { reversed: true, .. } => "reversed",
            Recipe::Fixed { .. } => "fixed",
            Recipe::Pairs { .. } => "pairs",
            Recipe::Ladder { .. } => "ladder",
            Recipe::Profmark { .. } => "profmark",
        }
    }
}

const USAGE: &str = "build-set --in <label>";

#[derive(Parser)]
#[command(name = "build-set", about = "cut a query and a target source into a set")]
struct Cli {
    /// Which label of paths.toml to build. Omit to list them
    #[arg(long = "in", value_name = "label")]
    label: Option<String>,

    /// Take back what is already at the label's `out` first, so a changed
    /// recipe can be rebuilt in place
    #[arg(long)]
    rebuild: bool,
}

fn main() -> anyhow::Result<()> {
    let paths = paths::File::open(env!("CARGO_MANIFEST_DIR"))?;

    let cli = Cli::parse();

    let Some(label) = cli.label.clone() else {
        println!("{}", listing(&paths)?);
        return Ok(());
    };

    let recipe: Recipe = paths.get(&label)?;

    if cli.rebuild {
        take_back(paths.at(recipe.out()))?;
    }

    match recipe {
        Recipe::Fixed {
            queries,
            alignments,
            targets,
            out,
            shards,
            seqs,
            fams,
            reversed,
            seed,
        } => fixed(
            Sources::new(&paths, queries, alignments, targets)?,
            paths.at(out),
            shards,
            seqs,
            fams,
            reversed,
            seed,
        ),
        Recipe::Pairs {
            queries,
            targets,
            out,
            pairs,
        } => make_pairs(
            paths.at(queries),
            paths.at(targets),
            paths.at(out),
            pairs,
        ),
        Recipe::Profmark {
            alignments,
            decoys,
            split,
            out,
            pairs,
            seed,
            train_test_id,
            min_test,
            max_test,
            threads,
        } => make_profmark(
            paths.at(alignments),
            paths.at(decoys),
            paths.at(split),
            paths.at(out),
            pairs,
            seed,
            train_test_id,
            min_test,
            max_test,
            threads,
        ),
        Recipe::Ladder {
            queries,
            alignments,
            targets,
            out,
            query_rungs,
            target_rungs,
            seed,
        } => ladder(
            Sources::new(&paths, queries, alignments, targets)?,
            paths.at(out),
            &query_rungs,
            &target_rungs,
            seed,
        ),
    }
}

/// The labels and what each would make, so the next command can be typed
/// without opening the file.
fn listing(paths: &paths::File) -> anyhow::Result<String> {
    use std::fmt::Write as _;

    let mut out = format!("labels in {}\n\n", paths.path().display());

    // the longest of each sets its column, so the three line up however the
    // labels are named and whatever a shape is called
    let width = paths.labels().map(str::len).max().unwrap_or(0);
    let shape_width = paths
        .labels()
        .map(|l| paths.get::<Recipe>(l).map(|r| r.shape().len()).unwrap_or(0))
        .max()
        .unwrap_or(0);

    for label in paths.labels() {
        let recipe: Recipe = paths.get(label)?;
        let (shape, size, dst) = match &recipe {
            Recipe::Fixed {
                out, shards, seqs, fams, ..
            } => (
                recipe.shape(),
                format!(
                    "{shards} {}, {} seqs, {} families",
                    plural(*shards, "shard"),
                    count(*seqs),
                    count(*fams)
                ),
                out,
            ),
            Recipe::Pairs { out, pairs, .. } => (
                recipe.shape(),
                match pairs {
                    Some(n) => format!("{n} {}", plural(*n, "pair")),
                    None => "every pair".to_string(),
                },
                out,
            ),
            Recipe::Profmark { out, pairs, .. } => (
                recipe.shape(),
                match pairs {
                    Some(n) => format!("{n} {}", plural(*n, "pair")),
                    None => "every pair that survives".to_string(),
                },
                out,
            ),
            Recipe::Ladder {
                out,
                query_rungs,
                target_rungs,
                ..
            } => (
                recipe.shape(),
                format!(
                    "{} query {} x {} target {}",
                    query_rungs.len(),
                    plural(query_rungs.len(), "rung"),
                    target_rungs.len(),
                    plural(target_rungs.len(), "rung")
                ),
                out,
            ),
        };

        let _ = writeln!(
            out,
            "  {label:<width$}  {shape:<shape_width$}  {size:<34}  -> {}",
            dst.display()
        );
    }

    let _ = write!(out, "\nusage: {USAGE}");
    Ok(out)
}

fn count(n: Option<usize>) -> String {
    match n {
        Some(n) => n.to_string(),
        None => "all".to_string(),
    }
}

fn shape_of(reversed: bool) -> &'static str {
    match reversed {
        true => set::shape::REVERSED.name,
        false => set::shape::FIXED.name,
    }
}

fn plural(n: usize, word: &str) -> String {
    match n {
        1 => word.to_string(),
        _ => format!("{word}s"),
    }
}

// -------------------------------------------------------------------- fixed

fn fixed(
    src: Sources,
    root: PathBuf,
    shards: usize,
    n_seqs: Option<usize>,
    n_fams: Option<usize>,
    reversed: bool,
    seed: u64,
) -> anyhow::Result<()> {
    claim(&root)?;

    let queries = root.join(QUERIES);
    let targets = root.join(TARGETS);
    let (query_hmm, query_sto) = (queries.join(QUERY_HMM), queries.join(QUERY_STO));

    // the deal counts as it writes, and the manifest is written once the
    // pipeline has finished, so what it counted comes back out through here
    let counted: Arc<Mutex<Vec<Shard>>> = Arc::default();

    PipelineBuilder::new()
        .step(
            PCmd::new("mkdir")
                .name("dirs")
                .flag("-p")
                .path(&queries)
                .path(&targets),
        )
        .step(
            Step::from_closures([
                Closure::new("query", {
                    let src = src.clone();
                    move || subset_query(&src, n_fams, &query_hmm, &query_sto)
                }),
                Closure::new("target", {
                    let src = src.clone();
                    let targets = targets.clone();
                    let counted = Arc::clone(&counted);

                    move || {
                        let seqs = src.collection()?;

                        let total = seqs.len();
                        let n_seqs = match n_seqs {
                            None => total,
                            Some(n) if n > total => {
                                eprintln!(
                                    "warning: asked for {n} sequences but the collection holds {total}"
                                );
                                total
                            }
                            Some(n) => n,
                        };

                        let dealt = deal(&seqs, n_seqs, shards, seed, reversed, &targets)?;
                        *counted.lock().expect("the deal poisoned the count") = dealt;
                        Ok(())
                    }
                }),
            ])
            .name("draw"),
        )
        .step(profile_db(&src.mmseqs, &queries).name("profile db"))
        .stderr_dir(root.join("stderr"))
        .sink(Progress::new())
        .build()?
        .run()?;

    let dealt = counted.lock().expect("the deal poisoned the count");
    let rows: Vec<set::Row> = dealt
        .iter()
        .map(|s| {
            set::Row::new(s.shard.to_string(), format!("{TARGETS}/{}.fa", s.shard))
                .query_hmm(format!("{QUERIES}/{QUERY_HMM}"))
                .query_sto(format!("{QUERIES}/{QUERY_STO}"))
                .query_db(format!("{QUERIES}/{QUERY_DB}"))
                .attr("shard", s.shard)
                .attr("seqs", s.seqs)
                .attr("residues", s.residues)
                .attr("bytes", s.bytes)
        })
        .collect();

    Set::new(&root, rows)
        .says("shape", shape_of(reversed))
        .says("recipe", "fixed")
        .says("seed", seed)
        .says("query", src.hmm.display())
        .says("targets", src.dir.display())
        .save()?;

    println!("\nbuilt {}", root.display());
    Ok(())
}

/// What one shard of a fixed set came to.
struct Shard {
    shard: usize,
    seqs: usize,
    residues: u64,
    bytes: u64,
}

/// Deals `n_seqs` sequences into `shards` files, `<i>.fa` for `i` in `1..=shards`.
///
/// The draw is a permutation, so plain round robin is enough to leave a shard
/// with nothing of the collection's own order in it.
fn deal(
    seqs: &Aggregate<IndexedFasta>,
    n_seqs: usize,
    shards: usize,
    seed: u64,
    reversed: bool,
    out_dir: &Path,
) -> anyhow::Result<Vec<Shard>> {
    std::fs::create_dir_all(out_dir)
        .with_context(|| format!("failed to create {}", out_dir.display()))?;

    let mut writers = Vec::with_capacity(shards);
    for i in 1..=shards {
        let path = out_dir.join(format!("{i}.fa"));
        let file =
            File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
        writers.push(BufWriter::new(file));
    }

    let mut rng = StdRng::seed_from_u64(seed);
    let drawn = seqs.permute_with(&mut rng).take(n_seqs);

    // counted as they are written: reading a thousand shards back to find out
    // how big they are is the whole deal again
    let mut counted = vec![(0usize, 0u64); shards];

    for (i, mut rec) in drawn.iter().enumerate() {
        let at = i % shards;

        // reversed as it is written rather than in a pass afterwards, so the
        // draw is the same one either way and a reversed set costs what a
        // forward one costs
        if reversed {
            rec.reverse();
        }

        rec.write_to(&mut writers[at], DEFAULT_LINE_WIDTH)?;

        counted[at].0 += 1;
        counted[at].1 += rec.seq.len() as u64;
    }

    for mut w in writers {
        w.flush()?;
    }

    counted
        .into_iter()
        .enumerate()
        .map(|(at, (seqs, residues))| {
            let shard = at + 1;
            let path = out_dir.join(format!("{shard}.fa"));
            let bytes = std::fs::metadata(&path)
                .with_context(|| format!("failed to stat {}", path.display()))?
                .len();

            Ok(Shard {
                shard,
                seqs,
                residues,
                bytes,
            })
        })
        .collect()
}

// -------------------------------------------------------------------- pairs

/// Pairs the `i`th query file with the `i`th target file, in name order.
///
/// A pair is the unit on both sides, so there is no draw and no seed: the two
/// directories already say what goes with what, and this copies them in and
/// writes down what each came to.
fn make_pairs(
    queries: PathBuf,
    targets: PathBuf,
    root: PathBuf,
    take: Option<usize>,
) -> anyhow::Result<()> {
    claim(&root)?;

    let (q, t) = (fastas(&queries)?, fastas(&targets)?);

    ensure!(
        q.len() == t.len(),
        "{} holds {} fasta files and {} holds {}; a pair needs one of each",
        queries.display(),
        q.len(),
        targets.display(),
        t.len()
    );

    let take = take.unwrap_or(q.len());
    ensure!(
        take <= q.len(),
        "asked for {take} pairs but there are {}",
        q.len()
    );

    let (q_dir, t_dir) = (root.join(QUERIES), root.join(TARGETS));
    std::fs::create_dir_all(&q_dir)?;
    std::fs::create_dir_all(&t_dir)?;

    let mut rows = Vec::new();
    for (pair, (from_q, from_t)) in q.iter().zip(&t).take(take).enumerate() {
        let pair = pair + 1;
        let (to_q, to_t) = (
            q_dir.join(format!("{pair}.fa")),
            t_dir.join(format!("{pair}.fa")),
        );

        std::fs::copy(from_q, &to_q)
            .with_context(|| format!("failed to copy {}", from_q.display()))?;
        std::fs::copy(from_t, &to_t)
            .with_context(|| format!("failed to copy {}", from_t.display()))?;

        let (q_size, t_size) = (measure(&to_q)?, measure(&to_t)?);

        rows.push(
            set::Row::new(pair.to_string(), format!("{TARGETS}/{pair}.fa"))
                .query_fa(format!("{QUERIES}/{pair}.fa"))
                .attr("pair", pair)
                .attr("query_residues", q_size.1)
                .attr("residues", t_size.1)
                .attr("seqs", t_size.0)
                .attr("bytes", std::fs::metadata(&to_t)?.len()),
        );
    }

    Set::new(&root, rows)
        .says("shape", set::shape::PAIRS.name)
        .says("recipe", "pairs")
        .says("queries", queries.display())
        .says("targets", targets.display())
        .save()?;

    println!("\nbuilt {}", root.display());
    Ok(())
}

/// The fasta a source names: every one in a directory, sorted by name, or the
/// single file itself.
//
// a source is a collection, and whether it arrived as one file or a thousand
// is the downloader's business rather than the recipe's. MGnify comes as
// shards and Swissprot as one file, and both are drawn from the same way
//
// sniffed rather than matched on the extension: these are
// somebody else's downloads, and a shard that arrived as
// .seq or with no extension at all is still a shard.
// detect_path errors on anything it cannot identify, which
// here is the same answer as "not a fasta"
fn fastas(dir: &Path) -> anyhow::Result<Vec<PathBuf>> {
    if dir.is_file() {
        return match format::detect_path(dir) {
            Ok(Format::Fasta) => Ok(vec![dir.to_path_buf()]),
            _ => bail!("{} is not a fasta", dir.display()),
        };
    }

    let mut out: Vec<PathBuf> = std::fs::read_dir(dir)
        .with_context(|| format!("failed to read {}", dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.is_file())
        .filter(|p| matches!(format::detect_path(p), Ok(Format::Fasta)))
        .collect();

    // read_dir hands them over in whatever order the filesystem
    // holds, and a record's position in the collection has to mean
    // the same thing on every run for a seed to reproduce a draw
    out.sort();

    if out.is_empty() {
        bail!("no fasta files in {}", dir.display());
    }

    Ok(out)
}

/// How many records a fasta holds and how many residues.
fn measure(path: &Path) -> anyhow::Result<(usize, u64)> {
    let fa = libsail::seq::fasta::Fasta::open(path)
        .with_context(|| format!("failed to read {}", path.display()))?;

    let residues = fa.iter().map(|rec| rec.seq.len() as u64).sum();

    Ok((fa.len(), residues))
}

// ------------------------------------------------------------------- ladder

fn ladder(
    src: Sources,
    root: PathBuf,
    query_rungs: &[usize],
    target_rungs: &[usize],
    seed: u64,
) -> anyhow::Result<()> {
    claim(&root)?;

    let queries = root.join(QUERIES);
    let targets = root.join(TARGETS);

    // ---- queries

    // an index of the file and nothing more: counting families does not need
    // the models parsed
    let n_fams = IndexedHmm::open(&src.hmm)
        .with_context(|| format!("failed to index {}", src.hmm.display()))?
        .len();
    let query_ladder = rungs(query_rungs, n_fams);
    println!("pfam holds {n_fams} families; query rungs: {query_ladder:?}");

    let mut query_residues: Vec<u64> = Vec::new();
    let mut pl = PipelineBuilder::new();

    for &q in &query_ladder {
        let dir = queries.join(q.to_string());
        std::fs::create_dir_all(&dir)
            .with_context(|| format!("failed to create {}", dir.display()))?;

        let (query_hmm, query_sto) = (dir.join(QUERY_HMM), dir.join(QUERY_STO));

        println!("taking {q} families...");
        // a rung at or past the whole of Pfam is the whole of Pfam, which the
        // subset says by being asked for nothing
        subset_query(&src, (q < n_fams).then_some(q), &query_hmm, &query_sto)?;

        // LENG is the query axis of a search's matrix, and so the honest
        // measure of how much work a rung is
        let models = IndexedHmm::open(&query_hmm)
            .with_context(|| format!("failed to index {}", query_hmm.display()))?;

        query_residues.push(models.iter().map(|m| m.header.leng as u64).sum());

        pl = pl.step(profile_db(&src.mmseqs, &dir).name(format!("queryDB.q{q}")));
    }

    println!("building the mmseqs profile dbs...");
    pl.stderr_dir(root.join("stderr"))
        .sink(Progress::new())
        .build()?
        .run()?;

    // ---- targets

    let seqs = src.collection()?;

    let total = seqs.len();
    let target_ladder = rungs(target_rungs, total);
    println!(
        "the collection holds {total} sequences across {} files; target rungs: {target_ladder:?}",
        seqs.parts().len()
    );

    let drawn = deal_nested(&seqs, &target_ladder, seed, &targets)?;

    // ---- the manifest

    // every rung of one axis against every rung of the other: a caller that
    // wants one query rung's targets filters on the attribute rather than
    // rebuilding the product
    let mut rows = Vec::new();
    for (&q, &q_residues) in query_ladder.iter().zip(&query_residues) {
        for (&t, &(t_residues, t_bytes)) in target_ladder.iter().zip(&drawn) {
            rows.push(
                set::Row::new(format!("q{q}.t{t}"), format!("{TARGETS}/{t}.fa"))
                    .query_hmm(format!("{QUERIES}/{q}/{QUERY_HMM}"))
                    .query_sto(format!("{QUERIES}/{q}/{QUERY_STO}"))
                    .query_db(format!("{QUERIES}/{q}/{QUERY_DB}"))
                    .attr("query_rung", q)
                    .attr("target_rung", t)
                    .attr("query_residues", q_residues)
                    .attr("target_residues", t_residues)
                    .attr("target_bytes", t_bytes),
            );
        }
    }

    Set::new(&root, rows)
        .says("shape", set::shape::LADDER.name)
        .says("recipe", "ladder")
        .says("seed", seed)
        .says("query", src.hmm.display())
        .says("targets", src.dir.display())
        .save()?;

    println!("\nbuilt {}", root.display());
    Ok(())
}

/// Sort the rungs, drop the duplicates, and cap them at what there is. Asking
/// for more than exists is how you say "all of it", so the rung is renamed to
/// the real number rather than refused.
fn rungs(asked: &[usize], max: usize) -> Vec<usize> {
    let mut out: Vec<usize> = asked
        .iter()
        .map(|&n| n.min(max))
        .filter(|&n| n > 0)
        .collect();

    out.sort_unstable();
    out.dedup();
    out
}

/// Draw every target rung in one pass, each a prefix of the next.
///
/// Returns the residues and bytes that landed in each.
fn deal_nested(
    seqs: &Aggregate<IndexedFasta>,
    rungs: &[usize],
    seed: u64,
    out_dir: &Path,
) -> anyhow::Result<Vec<(u64, u64)>> {
    std::fs::create_dir_all(out_dir)
        .with_context(|| format!("failed to create {}", out_dir.display()))?;

    let mut writers = Vec::with_capacity(rungs.len());
    for &n in rungs {
        let path = out_dir.join(format!("{n}.fa"));
        let file =
            File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
        writers.push(BufWriter::new(file));
    }

    let mut counted = vec![(0u64, 0u64); rungs.len()];

    let largest = rungs.last().copied().unwrap_or(0);
    let mut rng = StdRng::seed_from_u64(seed);
    let drawn = seqs.permute_with(&mut rng).take(largest);

    // the rungs are ascending, so record i belongs to every rung past the first
    // one big enough to hold it, and that boundary only ever moves forward
    let mut first = 0usize;

    // one record's bytes, written once and copied into every rung that holds
    // it -- and its length is the bytes that rung grew by
    let mut bytes: Vec<u8> = Vec::new();

    // every draw is a seek into a collection of billions, so the top rung takes
    // a while and says nothing while it does
    let start = std::time::Instant::now();
    const TICK: usize = 10_000;

    for (i, rec) in drawn.iter().enumerate() {
        while first < rungs.len() && rungs[first] <= i {
            first += 1;
        }

        bytes.clear();
        rec.write_to(&mut bytes, DEFAULT_LINE_WIDTH)?;
        let residues = rec.seq.len() as u64;

        for (w, c) in writers[first..].iter_mut().zip(counted[first..].iter_mut()) {
            w.write_all(&bytes)?;
            c.0 += residues;
            c.1 += bytes.len() as u64;
        }

        let done = i + 1;
        if done % TICK == 0 || done == largest {
            let secs = start.elapsed().as_secs_f64();
            let rate = done as f64 / secs;
            let left = (largest - done) as f64 / rate;
            eprint!("\r  drew {done}/{largest} ({rate:.0}/s, {left:.0}s left)    ");
        }
    }

    eprintln!();

    for mut w in writers {
        w.flush()?;
    }

    Ok(counted)
}


// ------------------------------------------------------------------- shared

/// What a build reads, and the binary it needs to finish.
#[derive(Clone)]
struct Sources {
    dir: std::path::PathBuf,
    hmm: std::path::PathBuf,
    sto: std::path::PathBuf,
    mmseqs: std::path::PathBuf,
}

impl Sources {
    /// What a label named, resolved against the file that named it.
    ///
    /// Nothing here defaults: a build reads what it was pointed at, so this
    /// crate holds no path to any download and cutting something other than
    /// Pfam is editing a line rather than changing code.
    fn new(
        paths: &paths::File,
        hmm: PathBuf,
        sto: PathBuf,
        targets: PathBuf,
    ) -> anyhow::Result<Sources> {
        Ok(Sources {
            dir: paths.at(targets),
            hmm: paths.at(hmm),
            sto: paths.at(sto),
            // checked up here: drawing the targets takes long enough that
            // finding out about a missing mmseqs afterwards would be miserable
            mmseqs: mmseqs()?,
        })
    }

    fn collection(&self) -> anyhow::Result<Aggregate<IndexedFasta>> {
        collection_at(&self.dir)
    }
}

/// Every fasta in `dir`, indexed and addressed as one collection.
///
/// The index is kept beside each file as `<name>.saidx`, so the pass over
/// every byte of the collection happens on the first build and not on the ones
/// after it. The index stamps the source's length and modification time, and
/// `index::is_current` refuses one that no longer matches, so a re-downloaded
/// shard is rebuilt rather than read through a stale map.
fn collection_at(dir: &Path) -> anyhow::Result<Aggregate<IndexedFasta>> {
    let paths = fastas(dir)?;

    let mut parts = Vec::with_capacity(paths.len());
    for path in &paths {
        parts.push(
            indexed(path).with_context(|| format!("failed to index {}", path.display()))?,
        );
    }

    Ok(Aggregate::new(parts))
}

/// One fasta, through its index on disk: opened if there is a current one,
/// built and written down if not.
///
/// Opened rather than read: `Index::open` takes the header and leaves the
/// offsets in the file, answering one at a time out of a `pread`. A draw of a
/// few hundred sequences touches a few hundred offsets, and reading the table
/// for every record of a 37 GB collection to reach them is the cost this
/// avoids.
///
/// A failure to write is a warning rather than an error. The index is a cache,
/// and a read-only or full source directory should cost the next build its
/// scan rather than this one its draw.
fn indexed(path: &Path) -> anyhow::Result<IndexedFasta> {
    let at = index::path_for(path);

    if index::is_current(&at, path).unwrap_or(false)
        && let Ok(index) = Index::open(&at)
    {
        return Ok(IndexedFasta::with_index(path, index)?);
    }

    let index = Index::build(File::open(path)?, Format::Fasta)?;

    if let Err(e) = index.write(&at, path) {
        eprintln!("warning: could not write an index beside {}: {e}", path.display());
    }

    Ok(IndexedFasta::with_index(path, index)?)
}

/// Refuse to build over an input set that is already there.
///
/// Rebuilding in place would leave whatever the last build wrote alongside
/// whatever this one does, and the pipelines read a directory rather than a
/// manifest, so the mixture would be searched as if it were one set.
/// Assemble a profmark benchmark and stamp the set it produced.
///
/// One unit: every query against one target file. The truth is per pair rather
/// than per unit, so `truth` names the file that carries it.
#[allow(clippy::too_many_arguments)]
fn make_profmark(
    alignments: PathBuf,
    decoys: PathBuf,
    split: PathBuf,
    out: PathBuf,
    pairs: Option<usize>,
    seed: u64,
    train_test_id: f64,
    min_test: usize,
    max_test: usize,
    threads: usize,
) -> anyhow::Result<()> {
    claim(&out)?;

    profmark::build(
        &alignments,
        &decoys,
        &split,
        &out,
        pairs,
        seed,
        train_test_id,
        min_test,
        max_test,
        threads,
        false,
    )?;

    let rows = vec![
        set::Row::new("1", "target.fa")
            .query_hmm("queries/query.hmm")
            .query_sto("queries/query.sto")
            .query_fa("queries/query.fa")
            .attr("truth", "truth.tbl"),
    ];

    Set::new(&out, rows)
        .says("shape", set::shape::PROFMARK.name)
        .says("recipe", "profmark")
        .says("seed", seed)
        .says("alignments", alignments.display())
        .says("decoys", decoys.display())
        .says("split", split.display())
        .save()?;

    println!("\nbuilt {}", out.display());
    Ok(())
}

fn claim(set: &Path) -> anyhow::Result<()> {
    if set.exists() {
        bail!(
            "{} already exists; pass --rebuild to take it back first",
            set.display()
        );
    }
    Ok(())
}

/// Remove what is at a label's `out`, so a changed recipe can be rebuilt.
//
// through the same measure-show-ask that `store clean` uses rather than a
// bare remove_dir_all: a fixed set at a thousand shards is most of a terabyte,
// and a recipe typo that points `out` somewhere unintended should be read
// before it is acted on rather than after
fn take_back(set: PathBuf) -> anyhow::Result<()> {
    if !set.exists() {
        return Ok(());
    }

    let parent = set.parent().unwrap_or(&set).to_owned();
    util::clean::run(&parent, &[("set", set)])
}

/// The first `n` families of Pfam, as both an hmm file and its alignments.
///
/// `None` is all of Pfam, where there is nothing to pick out and a copy beats
/// reading the whole file a line at a time.
fn subset_query(
    src: &Sources,
    n: Option<usize>,
    query_hmm: &Path,
    query_sto: &Path,
) -> anyhow::Result<()> {
    let Some(n) = n else {
        std::fs::copy(&src.hmm, query_hmm)
            .with_context(|| format!("failed to copy {}", src.hmm.display()))?;
        std::fs::copy(&src.sto, query_sto)
            .with_context(|| format!("failed to copy {}", src.sto.display()))?;
        return Ok(());
    };

    let names: HashSet<String> = cut::subset_hmm(&src.hmm, n, query_hmm)?;

    let kept = cut::subset_sto(&src.sto, &names, query_sto)?;
    if kept != names.len() {
        bail!(
            "kept {kept} stockholm records but the hmm subset named {}; \
             pfam.sto and pfam.hmm may be out of sync",
            names.len()
        );
    }

    Ok(())
}

/// The mmseqs profile db, built from the stockholm alignments beside it.
///
/// msaDB is an intermediate on the way to the profiles and is thrown away, so
/// it lands in the same directory rather than anywhere a pipeline would look.
fn profile_db(mmseqs: &Path, dir: &Path) -> Step {
    let msa_db = dir.join("msaDB");
    let query_db = dir.join("queryDB");

    Step::serial([
        PCmd::new("mkdir")
            .name("dirs")
            .flag("-p")
            .path(&msa_db)
            .path(&query_db),
        PCmd::new(mmseqs)
            .name("convertmsa")
            .sub("convertmsa")
            .arg("--identifier-field", 0)
            .path(dir.join("query.sto"))
            .path(msa_db.join("msaDB")),
        PCmd::new(mmseqs)
            .name("msa2profile")
            .sub("msa2profile")
            .arg("--match-mode", 1)
            .path(msa_db.join("msaDB"))
            .path(query_db.join("queryDB")),
        PCmd::new("rm").name("cleanup").flag("-rf").path(&msa_db),
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    /// A collection whose records name themselves, so a deal can be checked
    /// against what went into it.
    fn collection(name: &str, files: usize, per_file: usize) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("mgy-deal-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        for f in 0..files {
            let mut body = String::new();
            for r in 0..per_file {
                let len = 5 + (f * 7 + r * 13) % 40;
                body.push_str(&format!(">f{f}r{r}\n{}\n", "A".repeat(len)));
            }
            std::fs::write(dir.join(format!("{f:02}.fa")), body).unwrap();
        }

        dir
    }

    /// Deal into a fresh directory and read back what landed in each shard.
    fn shards_of(dir: &Path, n: usize, shards: usize) -> (Vec<Vec<String>>, PathBuf) {
        let agg = collection_at(dir).unwrap();
        let out = dir.join(format!("out-{n}-{shards}"));
        deal(&agg, n, shards, 67779, false, &out).unwrap();

        let names = (1..=shards)
            .map(|i| {
                let text = std::fs::read_to_string(out.join(format!("{i}.fa"))).unwrap();
                text.lines()
                    .filter(|l| l.starts_with('>'))
                    .map(|l| l[1..].split_whitespace().next().unwrap().to_string())
                    .collect()
            })
            .collect();

        (names, out)
    }

    fn assert_balanced(shards: &[Vec<String>]) {
        let sizes: Vec<usize> = shards.iter().map(|s| s.len()).collect();
        let lo = sizes.iter().min().unwrap();
        let hi = sizes.iter().max().unwrap();
        assert!(
            hi - lo <= 1,
            "shard sizes differ by more than one: {sizes:?}"
        );
    }

    #[test]
    fn dealing_everything_places_each_record_once() {
        let dir = collection("all", 3, 100);
        let (shards, _) = shards_of(&dir, 300, 7);

        assert_balanced(&shards);

        let mut got: Vec<String> = shards.into_iter().flatten().collect();
        assert_eq!(got.len(), 300, "every record should land somewhere");
        got.sort();
        got.dedup();
        assert_eq!(got.len(), 300, "no record should land twice");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn dealing_a_subset_draws_distinct_records() {
        let dir = collection("subset", 3, 100);
        let (shards, _) = shards_of(&dir, 90, 7);

        assert_balanced(&shards);

        let mut got: Vec<String> = shards.into_iter().flatten().collect();
        assert_eq!(got.len(), 90);
        got.sort();
        got.dedup();
        assert_eq!(got.len(), 90, "a subset must not repeat a record");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn shards_stay_balanced_when_the_count_does_not_divide() {
        let dir = collection("ragged", 2, 50);

        // 100 records over 7 shards leaves a remainder, as does the subset
        for (n, k) in [(100usize, 7usize), (97, 7), (13, 5), (100, 3)] {
            let (shards, out) = shards_of(&dir, n, k);
            assert_balanced(&shards);
            let placed: usize = shards.iter().map(|s| s.len()).sum();
            assert_eq!(placed, n, "n={n} k={k}: wrong number placed");
            std::fs::remove_dir_all(&out).ok();
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_subset_is_spread_across_the_collection() {
        let dir = collection("spread", 4, 100);
        let (shards, _) = shards_of(&dir, 200, 5);

        // the draw should reach every source file, not just the first ones
        let got: Vec<String> = shards.into_iter().flatten().collect();
        for f in 0..4 {
            assert!(
                got.iter().any(|n| n.starts_with(&format!("f{f}r"))),
                "no records drawn from file {f}"
            );
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn dealing_is_reproducible_for_a_seed() {
        let dir = collection("seed", 2, 60);
        let (a, out_a) = shards_of(&dir, 60, 5);
        std::fs::remove_dir_all(&out_a).ok();
        let (b, _) = shards_of(&dir, 60, 5);
        assert_eq!(a, b);

        std::fs::remove_dir_all(&dir).ok();
    }

    /// Every rung a prefix of the next, which is what makes the grid a surface.
    #[test]
    fn ladder_rungs_are_nested() {
        let dir = collection("ladder", 2, 60);
        let agg = collection_at(&dir).unwrap();

        let out = dir.join("ladder-out");
        let rungs = [10usize, 30, 90];
        let counted = deal_nested(&agg, &rungs, 67779, &out).unwrap();

        let read = |n: usize| std::fs::read(out.join(format!("{n}.fa"))).unwrap();
        let (small, mid, big) = (read(10), read(30), read(90));

        assert_eq!(&big[..small.len()], &small[..], "10 is not a prefix of 90");
        assert_eq!(&big[..mid.len()], &mid[..], "30 is not a prefix of 90");

        // and the counts describe what actually landed
        for (&rung, (_, bytes)) in rungs.iter().zip(&counted) {
            assert_eq!(*bytes, read(rung).len() as u64, "rung {rung}: wrong bytes");
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    /// Rungs past what exists collapse onto it rather than being refused.
    #[test]
    fn rungs_are_sorted_deduped_and_capped() {
        assert_eq!(rungs(&[100, 10, 100], 1000), vec![10, 100]);
        assert_eq!(rungs(&[10, 5000], 1000), vec![10, 1000]);
        assert_eq!(rungs(&[5000, 6000], 1000), vec![1000]);
        assert_eq!(rungs(&[0, 10], 1000), vec![10]);
    }
}
