use std::path::PathBuf;

use anyhow::{bail, Result};
use clap::Parser;

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

    let src_fa = repo.join("data/mgnify.fa");
    let src_hmm = repo.join("data/pfam.hmm");
    let src_sto = repo.join("data/pfam.sto");

    for path in [&src_fa, &src_hmm, &src_sto] {
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
    let source = match args.n_seqs {
        None => {
            println!("splitting all of MGnify into {} shards...", args.shards);
            src_fa.clone()
        }
        Some(n) => {
            println!("sampling {n} MGnify sequences...");
            let sampled = bench.join("target.fa");
            fasta::sample_to(&src_fa, n, &sampled)?;
            sampled
        }
    };

    println!("splitting into {} shards...", args.shards);
    fasta::split(&source, args.shards, &mgy, args.seed)?;

    // the sampled intermediate is redundant once it has been sharded
    if args.n_seqs.is_some() {
        std::fs::remove_file(bench.join("target.fa")).ok();
    }

    println!("\nbuilt {}", bench.display());
    Ok(())
}
