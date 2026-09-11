//! A recall directory the size of a real one, without the search.
//!
//! `mgy parse scores` is the thing being measured, and what it costs is set by
//! how many bytes of hit table it reads rather than by what produced them. So
//! this writes tables of the right shape and size -- the layouts the four
//! readers expect, hmmer's as a `cat` of query-split parts with its headers
//! repeated in the middle, MGnify names disjoint per shard -- and a
//! `ledger.tbl` over them.
//!
//! ```
//! synth_results --out outputs/synth --shards 10 --pairs 10000000
//! ```
//!
//! The pair sets are nested inside a tool: what a more sensitive run finds is
//! everything the less sensitive one found and more, and every run of one tool
//! gives a pair the same score. That is what the table's one column per tool
//! assumes, so a benchmark that broke it would not be measuring the same work.
//!
//! This is a dev tool. Nothing in a pipeline runs it.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::Context;
use clap::Parser;
use rayon::prelude::*;

use util::ledger::{Ledger, Row};

/// The runs, in the order the ledger declares them, with the fraction of a
/// shard's pairs each reports. Nested inside a tool.
const RUNS: [(&str, &str, f64); 6] = [
    ("nail-s9.0", "nail", 0.40),
    ("nail-s10.0", "nail", 0.55),
    ("nail-s12.0", "nail", 0.65),
    ("mmseqs-s7.5-ms2000", "mmseqs", 0.30),
    ("mmseqs-s12.0-ms2000", "mmseqs", 0.42),
    ("hmmer", "hmmer", 0.50),
];

/// How many rows one query-split part of an hmmer table holds, after which its
/// header turns up again.
const PART: usize = 50_000;

#[derive(Parser)]
#[command(about = "Write a recall directory of synthetic results")]
struct Args {
    /// Where to write the pipeline directory
    #[arg(long, value_name = "dir")]
    out: PathBuf,

    /// How many target shards
    #[arg(long, default_value_t = 10, value_name = "N")]
    shards: usize,

    /// How many query/target pairs per shard
    #[arg(long, default_value_t = 10_000_000, value_name = "N")]
    pairs: usize,

    /// How many families to draw query names from
    #[arg(long, default_value_t = 5000, value_name = "N")]
    fams: usize,

    /// The models to take family names from
    #[arg(long, value_name = "pfam.hmm")]
    pfam: Option<PathBuf>,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let names = families(args.pfam.as_deref(), args.fams)?;
    let results = args.out.join("results");
    std::fs::create_dir_all(&results)?;

    queries(&args.out.join("query.hmm"), &names)?;

    let shards: Vec<String> = (1..=args.shards).map(|n| n.to_string()).collect();
    let bytes: u64 = shards
        .par_iter()
        .map(|shard| shard_results(&results, shard, &names, args.pairs))
        .sum::<anyhow::Result<u64>>()?;

    targets(&args.out.join("targets"), &shards, args.pairs)?;
    ledger(&args.out, &shards)?;

    println!(
        "wrote {} ({:.1} GB over {} shards)",
        args.out.display(),
        bytes as f64 / (1u64 << 30) as f64,
        args.shards
    );

    Ok(())
}

/// Family names, read off a real `pfam.hmm` so the cutoffs table has something
/// to say about them.
fn families(pfam: Option<&Path>, want: usize) -> anyhow::Result<Vec<String>> {
    let path = match pfam {
        Some(path) => path.to_path_buf(),
        None => util::tools::pfam_hmm()?,
    };

    let file = File::open(&path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut names = Vec::with_capacity(want);

    use std::io::BufRead;
    for line in std::io::BufReader::with_capacity(1 << 20, file).lines() {
        let line = line?;
        if let Some(name) = line.strip_prefix("NAME  ") {
            names.push(name.trim().to_string());
            if names.len() == want {
                break;
            }
        }
    }

    names.sort_unstable();
    names.dedup();

    Ok(names)
}

/// A query set holding those families, so `parse` has one to number them from.
fn queries(path: &Path, names: &[String]) -> anyhow::Result<()> {
    let mut out = BufWriter::new(File::create(path)?);

    for (at, name) in names.iter().enumerate() {
        // a length that varies the way a real set's does
        out.write_all(util::profile(name, 40 + at % 300).as_bytes())?;
    }

    out.flush()?;
    Ok(())
}

/// What the shards would have come to, since nothing here writes the sequences
/// themselves.
fn targets(dir: &Path, shards: &[String], pairs: usize) -> anyhow::Result<()> {
    std::fs::create_dir_all(dir)?;

    let count = pairs / 10;
    let rows: String = shards
        .iter()
        .map(|shard| {
            format!(
                "  {shard} {count} {} {}\n",
                count as u64 * 320,
                count as u64 * 340
            )
        })
        .collect();

    std::fs::write(
        dir.join("sizes.tbl"),
        format!("# shard count residues bytes\n# ----- ----- -------- -----\n{rows}"),
    )?;

    Ok(())
}

fn ledger(dir: &Path, shards: &[String]) -> anyhow::Result<()> {
    let rows: Vec<Row> = RUNS
        .iter()
        .flat_map(|(name, tool, _)| {
            shards.iter().map(move |shard| Row {
                name: name.to_string(),
                tool: tool.to_string(),
                shard: shard.clone(),
                stage: String::new(),
                params: params(name),
                wall_s: Some(1.0),
            })
        })
        .collect();

    Ledger::from_rows(rows).write(&dir.join("ledger.tbl"))
}

/// The sensitivity out of a run's name, which is what tells two runs of one
/// tool apart in the summary.
fn params(name: &str) -> BTreeMap<String, String> {
    match name.split_once("-s") {
        Some((_, rest)) => BTreeMap::from([(
            "s".to_string(),
            rest.split('-').next().unwrap_or(rest).to_string(),
        )]),
        None => BTreeMap::new(),
    }
}

/// Every run's table for one shard.
fn shard_results(
    results: &Path,
    shard: &str,
    names: &[String],
    pairs: usize,
) -> anyhow::Result<u64> {
    // targets are dealt to one shard each, so a shard's names are its own
    let base = shard.parse::<u64>().unwrap_or(1) * 10_000_000_000;
    let targets = (pairs / 10).max(1) as u64;

    let mut rng = Rng::new(0x5eed ^ base);
    let mut written = 0u64;

    let mut tables: Vec<BufWriter<File>> = RUNS
        .iter()
        .map(|(name, ..)| {
            let path = results.join(format!("{name}.{shard}.tbl"));
            Ok(BufWriter::with_capacity(1 << 22, File::create(path)?))
        })
        .collect::<anyhow::Result<_>>()?;

    let mut dom = BufWriter::with_capacity(
        1 << 22,
        File::create(results.join(format!("hmmer.{shard}.domtbl")))?,
    );

    let hmmer = RUNS.len() - 1;
    let mut rows = [0usize; RUNS.len()];
    let mut dom_rows = 0usize;

    for _ in 0..pairs {
        let at = rng.next() % names.len() as u64;
        let number = base + rng.next() % targets;

        let query = &names[at as usize];
        let target = format!("MGYP{number:012}");

        // the pair decides its own scores, so drawing it twice writes the same
        // rows twice -- which is what a search that reported it twice would
        // have done, and what the table's one column per tool rests on
        let mut pair = Rng::new(at.wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ number);

        // how hard the pair is to find: a run reports it when its reach is
        // longer than this, which is what nests the sets inside a tool
        let hard = pair.unit();

        let nail = 5.0 + pair.unit() as f32 * 295.0;
        let mmseqs = nail * 0.9;
        let hmmer_score = nail * 1.02;

        for (at, (_, tool, reach)) in RUNS.iter().enumerate() {
            if hard > *reach {
                continue;
            }

            let head = header(&mut rows[at], *tool, PART);
            let table = &mut tables[at];

            match *tool {
                "nail" => {
                    table.write_all(head.as_bytes())?;
                    write!(
                        table,
                        "{target} {query} 1      100    1     100   {nail:.1}  0.0  1.2e-10 0.066\n"
                    )?;
                }
                "mmseqs" => {
                    table.write_all(head.as_bytes())?;
                    write!(
                        table,
                        "{query} {target} 31.0 100 60 2 1 100 1 100 1.2e-10 {mmseqs:.1}\n"
                    )?;
                }
                _ => {
                    let domains = 1 + (pair.next() % 3) as usize;
                    table.write_all(head.as_bytes())?;
                    write!(
                        table,
                        "{target}     -          {query:20} PF00000.1    1.2e-10 {hmmer_score:.1}   0.5   1.2e-10 {hmmer_score:.1}   0.5   1.1   1   0   0   1   {domains}   {domains}   {domains} FL=0\n"
                    )?;

                    for i in 0..domains {
                        let part = hmmer_score / domains as f32 - i as f32;
                        dom.write_all(header(&mut dom_rows, "domtbl", PART).as_bytes())?;
                        write!(
                            dom,
                            "{target}     -            178 {query:20} PF00000.1    232   1.2e-10 {hmmer_score:.1}   0.5   {} {domains}   1.1e-46   1.1e-43  {part:.1}   0.5     5   143    42   178    39   178 0.97 FL=0\n",
                            i + 1
                        )?;
                    }
                }
            }
        }

        let _ = hmmer;
    }

    for mut table in tables {
        table.flush()?;
        written += table.get_ref().metadata()?.len();
    }

    dom.flush()?;
    written += dom.get_ref().metadata()?.len();

    Ok(written)
}

/// hmmer's header, back again at the top of every part a `cat` joined.
fn header(rows: &mut usize, tool: &str, part: usize) -> &'static str {
    let at = *rows;
    *rows += 1;

    if tool == "nail" || tool == "mmseqs" || at % part != 0 {
        return "";
    }

    "#                                                               --- full sequence ---- --- best 1 domain ----\n\
     # target name        accession  query name           accession    E-value  score  bias\n\
     #------------------- ---------- -------------------- ---------- --------- ------ -----\n"
}

/// Enough randomness to spread the pairs out, and none of it worth a
/// dependency.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Rng {
        Rng(seed | 1)
    }

    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}
