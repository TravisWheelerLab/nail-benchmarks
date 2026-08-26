//! Draws the one query set and the target shards every pipeline searches.
//!
//! There is no benchmark-directory axis here the way there was in the crates
//! this replaces: the inputs live at the crate root and all three pipelines
//! read them, so a number from one is comparable to a number from another
//! without anyone having to check that the two draws happened to match.
//!
//! The stockholm alignments and the mmseqs profile db are always built, even
//! though only `recall` searches mmseqs today. They are what an mmseqs column
//! costs, and a query set that can't answer for one of the tools isn't one set.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, bail};
use clap::Parser;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;

use bioio::aggregate::AggregateFasta;
use bioio::{fasta, hmm, stockholm};
use feisty::Permutation;
use pail::{Closure, Cmd, PipelineBuilder, Progress, Step};
use tools::{mgnify, mmseqs, pfam_hmm, pfam_sto};

// the index that comes with the collection. building one from scratch is a pass
// over every byte of MGnify, so this is not a file to go without
const INDEX_NAME: &str = "mgnify.afi";

#[derive(Parser, Debug)]
pub struct Args {
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

pub fn main(args: Args) -> anyhow::Result<()> {
    let src_dir = mgnify()?;
    let src_hmm = pfam_hmm()?;
    let src_sto = pfam_sto()?;

    // checked up here: dealing the shards takes long enough that finding out
    // about a missing mmseqs afterwards would be miserable
    let mmseqs = mmseqs()?;

    let queries = crate::queries();
    let targets = crate::targets();

    for dir in [&queries, &targets] {
        if dir.exists() {
            bail!("{} already exists; remove it to rebuild", dir.display());
        }
    }

    let query_hmm = queries.join("query.hmm");
    let query_sto = queries.join("query.sto");

    let msa_db = queries.join("msaDB");
    let query_db = queries.join("queryDB");

    PipelineBuilder::new()
        .step(
            Cmd::new("mkdir")
                .name("dirs")
                .flag("-p")
                .path(&queries)
                .path(&targets),
        )
        .step(
            Step::from_closures([
                Closure::new("query", {
                    let (src_hmm, src_sto) = (src_hmm.clone(), src_sto.clone());
                    let (query_hmm, query_sto) = (query_hmm.clone(), query_sto.clone());
                    let n_fams = args.n_fams;

                    move || {
                        match n_fams {
                            // all of Pfam, so there is nothing to pick out, and
                            // a copy beats reading 1.6GB a line at a time
                            None => {
                                std::fs::copy(&src_hmm, &query_hmm).with_context(|| {
                                    format!("failed to copy {}", src_hmm.display())
                                })?;
                                std::fs::copy(&src_sto, &query_sto).with_context(|| {
                                    format!("failed to copy {}", src_sto.display())
                                })?;
                            }
                            Some(n) => {
                                let names = hmm::subset(&src_hmm, n, &query_hmm)?;

                                let kept = stockholm::subset_by_id(&src_sto, &names, &query_sto)?;
                                if kept != names.len() {
                                    bail!(
                                        "kept {kept} stockholm records but the hmm subset named \
                                         {}; pfam.sto and pfam.hmm may be out of sync",
                                        names.len()
                                    );
                                }
                            }
                        }

                        Ok(())
                    }
                }),
                Closure::new("target", {
                    let src_dir = src_dir.clone();
                    let (n_seqs, shards, seed, targets) =
                        (args.n_seqs, args.shards, args.seed, targets.clone());

                    move || {
                        // allow_overwrite only bites when the index no longer
                        // matches its sources, and then it rebuilds rather than
                        // refusing. that is a pass over every byte of MGnify,
                        // so a run that suddenly goes quiet for hours is this
                        let seqs = AggregateFasta::builder()
                            .dir(&src_dir)
                            .index(src_dir.join(INDEX_NAME))
                            .allow_overwrite()
                            .build()?;

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
        .step(
            Step::serial([
                Cmd::new("mkdir")
                    .name("dirs")
                    .flag("-p")
                    .path(&msa_db)
                    .path(&query_db),
                Cmd::new(&mmseqs)
                    .name("convertmsa")
                    .sub("convertmsa")
                    .arg("--identifier-field", 0)
                    .path(&query_sto)
                    .path(msa_db.join("msaDB")),
                Cmd::new(&mmseqs)
                    .name("msa2profile")
                    .sub("msa2profile")
                    .arg("--match-mode", 1)
                    .path(msa_db.join("msaDB"))
                    .path(query_db.join("queryDB")),
                Cmd::new("rm").name("cleanup").flag("-rf").path(&msa_db),
            ])
            .name("profile db"),
        )
        .stderr_dir(crate::dir().join("tmp/stderr"))
        .sink(Progress::new())
        .build()?
        .run()?;

    println!("\nbuilt {} and {}", queries.display(), targets.display());
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
}
