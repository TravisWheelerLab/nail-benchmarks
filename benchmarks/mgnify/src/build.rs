use std::collections::HashSet;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::Parser;

use crate::fasta;

/// This benchmark's directory, relative to the repository root.
pub const DIR: &str = "benchmarks/mgnify";

#[derive(Parser, Debug)]
pub struct Args {
    /// Names the output directory benchmark-<size>/.
    #[arg(short, long, default_value = "toy")]
    pub size: String,

    /// Number of target shards.
    #[arg(long, default_value_t = 4)]
    pub shards: usize,

    /// Target sequences to draw from MGnify; 0 uses all 2.4M.
    #[arg(long, default_value_t = 20_000)]
    pub seqs: usize,

    /// Pfam families to use as queries; 0 uses all of Pfam.
    #[arg(long, default_value_t = 50)]
    pub families: usize,

    /// Also write reversed copies of each shard, used as decoys when
    /// calibrating score cutoffs.
    #[arg(long)]
    pub reverse: bool,

    /// Seed for the shard deal.
    #[arg(long, default_value_t = 67779)]
    pub seed: u64,

    /// Repository root, if it cannot be discovered automatically.
    #[arg(long)]
    pub root: Option<PathBuf>,
}

pub fn main(args: Args) -> Result<()> {
    let repo = run::repo_root(args.root.as_deref())?;

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

    // everything the benchmark needs lives inside its own directory, so a size
    // can be rebuilt or deleted without touching data/ or any other size
    let bench = repo.join(DIR).join(format!("benchmark-{}", args.size));
    if bench.exists() {
        std::fs::remove_dir_all(&bench)?;
    }
    std::fs::create_dir_all(&bench)?;

    // ---- queries ----

    let query_hmm = bench.join("query.hmm");
    let query_sto = bench.join("query.sto");

    if args.families == 0 {
        println!("copying all of Pfam...");
        std::fs::copy(&src_hmm, &query_hmm)?;
        std::fs::copy(&src_sto, &query_sto)?;
    } else {
        println!("taking {} Pfam families...", args.families);
        let names = subset_hmm(&src_hmm, args.families, &query_hmm)?;
        if names.len() < args.families {
            bail!(
                "asked for {} families but Pfam yielded only {}",
                args.families,
                names.len()
            );
        }
        // the stockholm subset is driven by the names the hmm subset kept, so
        // the two always describe the same families
        subset_sto(&src_sto, &names, &query_sto)?;
    }

    // ---- targets ----

    let mgy = bench.join("mgy");
    let source = if args.seqs == 0 {
        println!("splitting all of MGnify into {} shards...", args.shards);
        src_fa.clone()
    } else {
        println!("sampling {} MGnify sequences...", args.seqs);
        let sampled = bench.join("target.fa");
        fasta::sample_to(&src_fa, args.seqs, &sampled)?;
        sampled
    };

    println!("splitting into {} shards...", args.shards);
    fasta::split(&source, args.shards, &mgy, args.seed)?;

    // the sampled intermediate is redundant once it has been sharded
    if args.seqs != 0 {
        std::fs::remove_file(bench.join("target.fa")).ok();
    }

    if args.reverse {
        println!("reversing {} shards...", args.shards);
        let rev = bench.join("mgy-rev");
        std::fs::create_dir_all(&rev)?;
        for i in 1..=args.shards {
            fasta::reverse(&mgy.join(format!("{i}.fa")), &rev.join(format!("{i}.fa")))
                .with_context(|| format!("failed to reverse shard {i}"))?;
        }
    }

    println!("\nbuilt {}", bench.display());
    Ok(())
}

/// Copy the first `n` HMM blocks of `src` into `dst`, returning their names.
///
/// Streams rather than parsing: Pfam is 1.6GB and only the block boundaries and
/// NAME lines matter here.
fn subset_hmm(src: &Path, n: usize, dst: &Path) -> Result<HashSet<String>> {
    use std::io::{BufRead, BufReader, BufWriter};

    let mut reader = BufReader::new(
        std::fs::File::open(src).with_context(|| format!("failed to open {}", src.display()))?,
    );
    let mut writer = BufWriter::new(std::fs::File::create(dst)?);

    let mut names = HashSet::new();
    let mut complete = 0usize;
    let mut line = String::new();

    // count completed blocks rather than names seen: stopping at the nth NAME
    // would truncate that block before its model and terminating //
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }

        if let Some(rest) = line.strip_prefix("NAME") {
            names.insert(rest.trim().to_string());
        }

        writer.write_all(line.as_bytes())?;

        if line.starts_with("//") {
            complete += 1;
            if complete == n {
                break;
            }
        }
    }

    writer.flush()?;
    Ok(names)
}

/// Copy the stockholm records whose `#=GF ID` is in `names` into `dst`.
fn subset_sto(src: &Path, names: &HashSet<String>, dst: &Path) -> Result<()> {
    use std::io::{BufRead, BufReader, BufWriter};

    let mut reader = BufReader::new(
        std::fs::File::open(src).with_context(|| format!("failed to open {}", src.display()))?,
    );
    let mut writer = BufWriter::new(std::fs::File::create(dst)?);

    let mut block: Vec<String> = Vec::new();
    let mut id: Option<String> = None;
    let mut kept = 0usize;
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }

        if let Some(rest) = line.strip_prefix("#=GF ID") {
            id = Some(rest.trim().to_string());
        }

        block.push(line.clone());

        if line.starts_with("//") {
            if id.as_ref().is_some_and(|i| names.contains(i)) {
                for l in &block {
                    writer.write_all(l.as_bytes())?;
                }
                kept += 1;
                if kept == names.len() {
                    break;
                }
            }
            block.clear();
            id = None;
        }
    }

    writer.flush()?;

    if kept != names.len() {
        bail!(
            "kept {kept} stockholm records but the hmm subset named {}; \
             pfam.sto and pfam.hmm may be out of sync",
            names.len()
        );
    }

    Ok(())
}
