use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{bail, Context};
use clap::Parser;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;

use bioio::aggregate::AggregateFasta;
use bioio::{fasta, hmm, stockholm};
use feisty::Permutation;
use pipeline::{Cmd, PipelineBuilder, Progress, Step};
use tools::{mgnify, mmseqs, pfam_hmm, pfam_sto};

pub const DEFAULT_NAME: &str = "benchmark";
const INDEX_NAME: &str = "mgnify.afi";

#[derive(Parser, Debug)]
pub struct Args {
    /// The name of the benchmark directory. The path resolves to "benchmarks/mgnify/<name>/"
    #[arg(default_value = DEFAULT_NAME)]
    pub name: String,

    /// The number of target database shards
    #[arg(long, default_value_t = 1000, value_name = "N")]
    pub shards: usize,

    /// Impose a limit on the number of MGnify sequences used for the benchmark
    #[arg(long = "seqs", value_name = "N")]
    pub n_seqs: Option<usize>,

    /// Impose a limit on the number of Pfam families used for the benchmark
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

    let mmseqs = mmseqs()?;

    let bench = crate::util::dir().join(&args.name);

    if bench.exists() {
        bail!("benchmark: {} already exists", args.name)
    }

    std::fs::create_dir_all(&bench)?;

    // ---- queries ----

    let queries = bench.join("queries");
    std::fs::create_dir_all(&queries)?;

    let query_hmm = queries.join("query.hmm");
    let query_sto = queries.join("query.sto");

    match args.n_fams {
        None => {
            println!("copying all of Pfam...");
            std::fs::copy(&src_hmm, &query_hmm)?;
            std::fs::copy(&src_sto, &query_sto)?;
        }
        Some(n) => {
            println!("taking {n} Pfam families...");
            let names = hmm::subset(&src_hmm, n, &query_hmm)?;

            let kept = stockholm::subset_by_id(&src_sto, &names, &query_sto)?;
            if kept != names.len() {
                bail!(
                    "kept {kept} stockholm records but the hmm subset named {}; \
                     pfam.sto and pfam.hmm may be out of sync",
                    names.len()
                );
            }
        }
    }

    println!("building the mmseqs profile db...");

    let msa_db = queries.join("msaDB");
    let query_db = queries.join("queryDB");

    PipelineBuilder::new()
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
        .stderr_dir(bench.join("stderr"))
        .sink(Progress::new())
        .build()
        .run()?;

    // ---- targets ----

    let targets = bench.join("targets");

    let seqs = AggregateFasta::builder()
        .dir(&src_dir)
        .index(src_dir.join(INDEX_NAME))
        .allow_overwrite()
        .build()?;

    let total = seqs.len();
    let n_seqs = match args.n_seqs {
        None => total,
        Some(n) if n as u64 > total => {
            eprintln!("warning: asked for {n} sequences but the collection holds {total}");
            total
        }
        Some(n) => n as u64,
    };

    println!(
        "dealing {n_seqs} of {total} sequences across {} files into {} shards...",
        seqs.files().len(),
        args.shards
    );
    deal(&seqs, n_seqs, args.shards, args.seed, &targets)?;

    println!("\nbuilt {}", bench.display());
    Ok(())
}

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
        let dir = std::env::temp_dir().join(format!("mgnify-deal-{}-{name}", std::process::id()));
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
