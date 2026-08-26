//! Cuts Pfam and MGnify into an input set for the pipelines to search.
//!
//! There is no benchmark-directory axis here the way there was in the crates
//! this replaces: an input set lives at the crate root under its kind, and
//! every pipeline of that kind reads it, so a number from one is comparable to
//! a number from another without anyone having to check that the two draws
//! happened to match.
//!
//! The two kinds differ in how they cut, not in what they cut -- see
//! [`crate::inputs`] for the shapes. What they share is here: subsetting the
//! query, and building the mmseqs profile db from what came out.
//!
//! The stockholm alignments and the mmseqs profile db are always built, even
//! for the pipelines that don't search mmseqs. They are what an mmseqs column
//! costs, and a query set that can't answer for one of the tools isn't one set.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, bail};
use clap::{Parser, Subcommand};
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;

use bioio::aggregate::AggregateFasta;
use bioio::split::{self, Kind as SplitKind};
use bioio::{fasta, hmm, stockholm};
use feisty::Permutation;
use pail::{Closure, Cmd as PCmd, PipelineBuilder, Progress, Step};
use tools::{mgnify, mmseqs, pfam_hmm, pfam_sto};

use crate::inputs;

// the index that comes with the collection. building one from scratch is a pass
// over every byte of MGnify, so this is not a file to go without
const INDEX_NAME: &str = "mgnify.afi";

#[derive(Subcommand)]
pub enum Cmd {
    /// One query set against target shards of equal size.
    Fixed(FixedArgs),
    /// Nested rungs on both axes, each a prefix of the next one up.
    Ladder(LadderArgs),
}

pub fn main(cmd: Cmd) -> anyhow::Result<()> {
    match cmd {
        Cmd::Fixed(args) => fixed(args),
        Cmd::Ladder(args) => ladder(args),
    }
}

// -------------------------------------------------------------------- fixed

#[derive(Parser, Debug)]
pub struct FixedArgs {
    /// The number of target shards
    #[arg(long, default_value_t = 1000, value_name = "N")]
    pub shards: usize,

    /// Impose a limit on the number of MGnify sequences used
    #[arg(long = "seqs", value_name = "N")]
    pub n_seqs: Option<usize>,

    /// Impose a limit on the number of Pfam families used
    #[arg(long = "fams", value_name = "N")]
    pub n_fams: Option<usize>,

    /// Random seed
    #[arg(long, default_value_t = 67779, value_name = "N")]
    pub seed: u64,
}

fn fixed(args: FixedArgs) -> anyhow::Result<()> {
    let src = Sources::find()?;

    let queries = inputs::fixed::queries();
    let targets = inputs::fixed::targets();
    claim(&[&queries, &targets])?;

    let (query_hmm, query_sto) = (inputs::fixed::query_hmm(), inputs::fixed::query_sto());

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
                    let (src, n_fams) = (src.clone(), args.n_fams);
                    move || subset_query(&src, n_fams, &query_hmm, &query_sto)
                }),
                Closure::new("target", {
                    let src = src.clone();
                    let (n_seqs, shards, seed) = (args.n_seqs, args.shards, args.seed);

                    move || {
                        let seqs = src.collection()?;

                        let total = seqs.len();
                        let n_seqs = match n_seqs {
                            None => total,
                            Some(n) if n as u64 > total => {
                                eprintln!(
                                    "warning: asked for {n} sequences but the collection holds {total}"
                                );
                                total
                            }
                            Some(n) => n as u64,
                        };

                        deal(&seqs, n_seqs, shards, seed, &targets)
                    }
                }),
            ])
            .name("draw"),
        )
        .step(profile_db(&src.mmseqs, &queries).name("profile db"))
        .stderr_dir(crate::dir().join("tmp/stderr"))
        .sink(Progress::new())
        .build()?
        .run()?;

    println!(
        "\nbuilt {} and {}",
        inputs::fixed::queries().display(),
        inputs::fixed::targets().display()
    );
    Ok(())
}

/// Deals `n_seqs` sequences into `shards` files, `<i>.fa` for `i` in `1..=shards`.
///
/// Dealing everything reshuffles the shard order every round rather than
/// drawing a permutation, since a permutation over the whole collection is a
/// cost this doesn't need to pay when nothing is actually being left out.
fn deal(
    seqs: &AggregateFasta,
    n_seqs: u64,
    shards: usize,
    seed: u64,
    out_dir: &Path,
) -> anyhow::Result<()> {
    std::fs::create_dir_all(out_dir)
        .with_context(|| format!("failed to create {}", out_dir.display()))?;

    let mut writers = Vec::with_capacity(shards);
    for i in 1..=shards {
        let path = out_dir.join(format!("{i}.fa"));
        let file =
            File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
        writers.push(BufWriter::new(file));
    }

    if n_seqs == seqs.len() {
        // sequences arrive in collection order, so the shard order is
        // reshuffled every round to break up neighbours
        let mut rng = StdRng::seed_from_u64(seed);
        let mut order: Vec<usize> = (0..shards).collect();

        let mut i = 0;
        for path in seqs.files() {
            let mut reader = fasta::Reader::from_path(path)?;
            while let Some(rec) = reader.next_record()? {
                let slot = i % shards;
                if slot == 0 {
                    order.shuffle(&mut rng);
                }
                writeln!(&mut writers[order[slot]], "{rec}")?;
                i += 1;
            }
        }
    } else {
        let perm = Permutation::new(seqs.len(), seed);
        let mut records = seqs.records();

        for i in 0..n_seqs {
            let bytes = records.get(perm.get(i))?;
            writers[(i % shards as u64) as usize].write_all(&bytes)?;
        }
    }

    for mut w in writers {
        w.flush()?;
    }

    Ok(())
}

// ------------------------------------------------------------------- ladder

#[derive(Parser, Debug)]
pub struct LadderArgs {
    /// Query rungs, as Pfam family counts. Anything at or past the number of
    /// families Pfam holds becomes all of Pfam.
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "10,100,1000,10000,100000",
        value_name = "N,N,..."
    )]
    pub queries: Vec<usize>,

    /// Target rungs, as MGnify sequence counts.
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "1000,10000,100000",
        value_name = "N,N,..."
    )]
    pub targets: Vec<usize>,

    /// Random seed
    #[arg(long, default_value_t = 67779, value_name = "N")]
    pub seed: u64,
}

fn ladder(args: LadderArgs) -> anyhow::Result<()> {
    let src = Sources::find()?;

    let queries = inputs::ladder::queries();
    let targets = inputs::ladder::targets();
    claim(&[&queries, &targets])?;

    // ---- queries

    let n_fams = split::index(&src.hmm, SplitKind::Hmm)?.len();
    let ladder = rungs(&args.queries, n_fams);
    println!("pfam holds {n_fams} families; query rungs: {ladder:?}");

    let mut sizes: Vec<Size> = Vec::new();
    let mut pl = PipelineBuilder::new();

    for &q in &ladder {
        let dir = inputs::ladder::query(q);
        std::fs::create_dir_all(&dir)
            .with_context(|| format!("failed to create {}", dir.display()))?;

        let (query_hmm, query_sto) = (inputs::ladder::query_hmm(q), inputs::ladder::query_sto(q));

        println!("taking {q} families...");
        // a rung at or past the whole of Pfam is the whole of Pfam, which the
        // subset says by being asked for nothing
        subset_query(&src, (q < n_fams).then_some(q), &query_hmm, &query_sto)?;

        // LENG per model, which is what the index weights records by
        let models = split::index(&query_hmm, SplitKind::Hmm)?;
        sizes.push(Size {
            rung: q,
            residues: models.iter().map(|r| r.weight).sum(),
            bytes: std::fs::metadata(&query_hmm)?.len(),
        });

        pl = pl.step(profile_db(&src.mmseqs, &dir).name(format!("queryDB.q{q}")));
    }

    write_sizes(&inputs::ladder::sizes(&queries), &sizes)?;

    println!("building the mmseqs profile dbs...");
    pl.stderr_dir(crate::dir().join("tmp/stderr"))
        .sink(Progress::new())
        .build()?
        .run()?;

    // ---- targets

    let seqs = src.collection()?;

    let total = seqs.len();
    let ladder = rungs(&args.targets, total as usize);
    println!(
        "the collection holds {total} sequences across {} files; target rungs: {ladder:?}",
        seqs.files().len()
    );

    let drawn = deal_nested(&seqs, &ladder, args.seed, &targets)?;

    let sizes: Vec<Size> = ladder
        .iter()
        .zip(drawn)
        .map(|(&rung, (residues, bytes))| Size {
            rung,
            residues,
            bytes,
        })
        .collect();

    write_sizes(&inputs::ladder::sizes(&targets), &sizes)?;

    println!("\nbuilt {} and {}", queries.display(), targets.display());
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
    seqs: &AggregateFasta,
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

    let largest = rungs.last().copied().unwrap_or(0) as u64;
    let perm = Permutation::new(seqs.len(), seed);
    let mut records = seqs.records();

    // the rungs are ascending, so record i belongs to every rung past the first
    // one big enough to hold it, and that boundary only ever moves forward
    let mut first = 0usize;

    // every draw is a seek into a collection of billions, so the top rung takes
    // a while and says nothing while it does
    let start = std::time::Instant::now();
    const TICK: u64 = 10_000;

    for i in 0..largest {
        while first < rungs.len() && (rungs[first] as u64) <= i {
            first += 1;
        }

        let bytes = records.get(perm.get(i))?;
        let residues = count_residues(&bytes);

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

/// Residues in one fasta record: everything past the header line that isn't
/// whitespace.
fn count_residues(record: &[u8]) -> u64 {
    let body = match record.iter().position(|&b| b == b'\n') {
        Some(i) => &record[i + 1..],
        None => return 0,
    };

    body.iter().filter(|b| !b.is_ascii_whitespace()).count() as u64
}

struct Size {
    rung: usize,
    residues: u64,
    bytes: u64,
}

/// What each rung actually turned out to be, since sequences and models are not
/// uniform units of work and residues is the honest axis to plot against.
///
/// Only the ladder writes this. The fixed kind's sizes travel in `scores.tbl`,
/// measured at parse time from the files that are present, and a second copy
/// here would be a second thing to keep true.
fn write_sizes(path: &Path, rows: &[Size]) -> anyhow::Result<()> {
    let headers = ["rung", "residues", "bytes"];

    let cells: Vec<[String; 3]> = rows
        .iter()
        .map(|r| {
            [
                r.rung.to_string(),
                r.residues.to_string(),
                r.bytes.to_string(),
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

    out.push('#');
    for (h, &w) in headers.iter().zip(&widths) {
        out.push_str(&format!(" {h:<w$}"));
    }
    out.push_str("\n#");
    for &w in &widths {
        out.push_str(&format!(" {}", "-".repeat(w)));
    }
    out.push('\n');

    for row in &cells {
        // the two the `# ` takes on a header line, so the columns sit under
        // their names rather than beside them
        out.push_str("  ");
        for (i, (c, &w)) in row.iter().zip(&widths).enumerate() {
            if i > 0 {
                out.push(' ');
            }
            out.push_str(&format!("{c:<w$}"));
        }
        out.push('\n');
    }

    std::fs::write(path, out).with_context(|| format!("failed to write {}", path.display()))
}

// ------------------------------------------------------------------- shared

/// The three files every build reads, and the binary it needs to finish.
#[derive(Clone)]
struct Sources {
    dir: std::path::PathBuf,
    hmm: std::path::PathBuf,
    sto: std::path::PathBuf,
    mmseqs: std::path::PathBuf,
}

impl Sources {
    fn find() -> anyhow::Result<Sources> {
        Ok(Sources {
            dir: mgnify()?,
            hmm: pfam_hmm()?,
            sto: pfam_sto()?,
            // checked up here: drawing the targets takes long enough that
            // finding out about a missing mmseqs afterwards would be miserable
            mmseqs: mmseqs()?,
        })
    }

    /// allow_overwrite only bites when the index no longer matches its sources,
    /// and then it rebuilds rather than refusing. That is a pass over every byte
    /// of MGnify, so a run that suddenly goes quiet for hours is this.
    fn collection(&self) -> anyhow::Result<AggregateFasta> {
        AggregateFasta::builder()
            .dir(&self.dir)
            .index(self.dir.join(INDEX_NAME))
            .allow_overwrite()
            .build()
    }
}

/// Refuse to build over an input set that is already there.
///
/// Rebuilding in place would leave whatever the last build wrote alongside
/// whatever this one does, and the pipelines read a directory rather than a
/// manifest, so the mixture would be searched as if it were one set.
fn claim(dirs: &[&Path]) -> anyhow::Result<()> {
    for dir in dirs {
        if dir.exists() {
            bail!("{} already exists; remove it to rebuild", dir.display());
        }
    }
    Ok(())
}

/// The first `n` families of Pfam, as both an hmm file and its alignments.
///
/// `None` is all of Pfam, where there is nothing to pick out and a copy beats
/// reading 1.6GB a line at a time.
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

    let names: HashSet<String> = hmm::subset(&src.hmm, n, query_hmm)?;

    let kept = stockholm::subset_by_id(&src.sto, &names, query_sto)?;
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
    fn shards_of(dir: &Path, n: u64, shards: usize) -> (Vec<Vec<String>>, PathBuf) {
        let agg = AggregateFasta::builder().dir(dir).build().unwrap();
        let out = dir.join(format!("out-{n}-{shards}"));
        deal(&agg, n, shards, 67779, &out).unwrap();

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
        for (n, k) in [(100u64, 7usize), (97, 7), (13, 5), (100, 3)] {
            let (shards, out) = shards_of(&dir, n, k);
            assert_balanced(&shards);
            let placed: usize = shards.iter().map(|s| s.len()).sum();
            assert_eq!(placed as u64, n, "n={n} k={k}: wrong number placed");
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
        let agg = AggregateFasta::builder().dir(&dir).build().unwrap();

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
