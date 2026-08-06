use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{bail, Result};
use clap::Parser;

use bioio::aggregate::AggregateFasta;
use bioio::{fasta, hmm, stockholm};

/// This benchmark's directory, fixed at compile time.
pub fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The repository root.
pub fn repo() -> PathBuf {
    run::repo(env!("CARGO_MANIFEST_DIR"))
}

/// Name of the benchmark directory when none is given.
pub const DEFAULT_NAME: &str = "benchmark";

/// Where the collection index is cached, beside the source fastas.
const INDEX_NAME: &str = "mgnify.afi";

#[derive(Parser, Debug)]
pub struct Args {
    /// Directory to build into, under benchmarks/mgnify/. Use this to keep
    /// several benchmarks side by side, e.g. a small one for development.
    #[arg(long, default_value = DEFAULT_NAME)]
    pub name: String,

    /// Number of target shards.
    #[arg(long, default_value_t = 1000)]
    pub shards: usize,

    /// Target sequences to draw from MGnify. Omit to use all of it; set it to
    /// build something small enough to develop against.
    #[arg(long = "seqs")]
    pub n_seqs: Option<usize>,

    /// Pfam families to use as queries. Omit to use all of Pfam.
    #[arg(long = "fams")]
    pub n_fams: Option<usize>,

    /// Seed for the shard deal.
    #[arg(long, default_value_t = 67779)]
    pub seed: u64,
}

pub fn main(args: Args) -> Result<()> {
    let repo = repo();

    let src_dir = repo.join("data/mgnify");
    let src_hmm = repo.join("data/pfam.hmm");
    let src_sto = repo.join("data/pfam.sto");

    for path in [&src_dir, &src_hmm, &src_sto] {
        if !path.exists() {
            bail!(
                "missing source data {}; run `make setup` from the repo root",
                path.display()
            );
        }
    }

    // everything for a benchmark — inputs, results, and analysis — lives under
    // one directory, so it can be rebuilt or deleted as a unit
    let bench = dir().join(&args.name);
    if bench.exists() {
        std::fs::remove_dir_all(&bench)?;
    }
    std::fs::create_dir_all(&bench)?;

    // ---- queries ----

    let query_hmm = bench.join("query.hmm");
    let query_sto = bench.join("query.sto");

    match args.n_fams {
        None => {
            println!("copying all of Pfam...");
            std::fs::copy(&src_hmm, &query_hmm)?;
            std::fs::copy(&src_sto, &query_sto)?;
        }
        Some(n) => {
            println!("taking {n} Pfam families...");
            let names = hmm::subset(&src_hmm, n, &query_hmm)?;

            // the stockholm subset is driven by the names the hmm subset kept, so
            // the two always describe the same families
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

    // ---- targets ----

    let mgy = bench.join("mgy");
    let sampled = bench.join("target.fa");

    match args.n_seqs {
        Some(n) => {
            // indexing the whole collection is the expensive part, so the index
            // is kept beside the source and reused by later builds
            let collection = AggregateFasta::builder()
                .dir(&src_dir)
                .index(src_dir.join(INDEX_NAME))
                .allow_overwrite()
                .build()?;

            println!(
                "sampling {n} of {} sequences across {} files...",
                collection.len(),
                collection.files().len()
            );

            let mut out = BufWriter::new(File::create(&sampled)?);
            collection.sample(n, args.seed, &mut out)?;
            out.flush()?;
        }
        None => todo!("sharding the whole collection needs the streaming splitter"),
    }

    println!("splitting into {} shards...", args.shards);
    fasta::split(&sampled, args.shards, &mgy, args.seed)?;

    // the sampled intermediate is redundant once it has been sharded
    std::fs::remove_file(&sampled).ok();

    println!("\nbuilt {}", bench.display());
    Ok(())
}
