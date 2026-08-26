//! Draws the two ladders a search-size run sweeps over.
//!
//! Both are nested: the query set at one rung is a prefix of the next one up,
//! and so is the target set. That makes the grid a surface rather than a pile of
//! unrelated points, and it means the whole target ladder comes out of one pass
//! over MGnify.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use anyhow::{bail, Context};
use clap::Parser;

use bioio::aggregate::AggregateFasta;
use bioio::split::{self, Kind};
use feisty::Permutation;
use pail::{Cmd, PipelineBuilder, Progress, Step};
use tools::{mgnify, mmseqs, pfam_hmm, pfam_sto};

pub const DEFAULT_NAME: &str = "benchmark";
const INDEX_NAME: &str = "mgnify.afi";

#[derive(Parser, Debug)]
pub struct Args {
    /// The name of the benchmark directory. The path resolves to "benchmarks/search-size/<name>/"
    #[arg(default_value = DEFAULT_NAME)]
    pub name: String,

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

pub fn main(args: Args) -> anyhow::Result<()> {
    let src_dir = mgnify()?;
    let src_hmm = pfam_hmm()?;
    let src_sto = pfam_sto()?;

    // checked up here: drawing the target ladder takes long enough that finding
    // out about a missing mmseqs afterwards would be miserable
    let mmseqs = mmseqs()?;

    let bench = crate::dir().join(&args.name);

    if bench.exists() {
        bail!("benchmark: {} already exists", args.name)
    }

    std::fs::create_dir_all(&bench)?;

    let mut sizes: Vec<Row> = Vec::new();

    // ---- queries ----

    let n_fams = count_models(&src_hmm)?;
    let queries = ladder(&args.queries, n_fams);
    println!("pfam holds {n_fams} families; query rungs: {queries:?}");

    let mut pl = PipelineBuilder::new();

    for &q in &queries {
        let dir = bench.join(format!("queries/{q}"));
        std::fs::create_dir_all(&dir)?;

        let query_hmm = dir.join("query.hmm");
        let query_sto = dir.join("query.sto");

        // the top rung is all of Pfam, so there is nothing to pick out
        if q == n_fams {
            println!("copying all of Pfam...");
            std::fs::copy(&src_hmm, &query_hmm)?;
            std::fs::copy(&src_sto, &query_sto)?;
        } else {
            println!("taking {q} families...");
            let names = subset_hmm(&src_hmm, q, &query_hmm)?;

            let kept = subset_sto(&src_sto, &names, &query_sto)?;
            if kept != names.len() {
                bail!(
                    "kept {kept} stockholm records but the hmm subset named {}; \
                     pfam.sto and pfam.hmm may be out of sync",
                    names.len()
                );
            }
        }

        // LENG per model, which is what the index weights records by
        let positions: u64 = split::index(&query_hmm, Kind::Hmm)?
            .iter()
            .map(|r| r.weight)
            .sum();

        sizes.push(Row {
            axis: "query",
            rung: q,
            residues: positions,
            bytes: std::fs::metadata(&query_hmm)?.len(),
        });

        let msa_db = dir.join("msaDB");
        let query_db = dir.join("queryDB");

        pl = pl.step(
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
            .name(format!("queryDB.q{q}")),
        );
    }

    println!("building the mmseqs profile dbs...");
    pl.stderr_dir(bench.join("stderr"))
        .sink(Progress::new())
        .build()?
        .run()?;

    // ---- targets ----

    let seqs = AggregateFasta::builder()
        .dir(&src_dir)
        .index(src_dir.join(INDEX_NAME))
        .allow_overwrite()
        .build()?;

    let total = seqs.len();
    let targets = ladder(&args.targets, total as usize);
    println!(
        "the collection holds {total} sequences across {} files; target rungs: {targets:?}",
        seqs.files().len()
    );

    let drawn = deal(&seqs, &targets, args.seed, &bench.join("targets"))?;
    for (&t, (residues, bytes)) in targets.iter().zip(drawn) {
        sizes.push(Row {
            axis: "target",
            rung: t,
            residues,
            bytes,
        });
    }

    write_sizes(&bench.join("sizes.tbl"), &sizes)?;

    println!("\nbuilt {}", bench.display());
    Ok(())
}

/// Sort the rungs, drop the duplicates, and cap them at what there is. Asking
/// for more than exists is how you say "all of it", so the rung is renamed to
/// the real number rather than refused.
fn ladder(rungs: &[usize], max: usize) -> Vec<usize> {
    let mut out: Vec<usize> = rungs
        .iter()
        .map(|&n| n.min(max))
        .filter(|&n| n > 0)
        .collect();

    out.sort_unstable();
    out.dedup();
    out
}

/// Copy the first `n` complete models of `src` into `dst`, returning their
/// names.
///
/// Completion is counted by `//`, not by `NAME`: stopping at the nth name would
/// truncate that model before its emission lines and terminator.
///
/// This works on bytes rather than lines because Pfam isn't valid utf-8 the
/// whole way through, and a line reader gives up when it reaches that. Names
/// themselves are ascii, so reading them lossily costs nothing.
fn subset_hmm(src: &Path, n: usize, dst: &Path) -> anyhow::Result<HashSet<String>> {
    let file = File::open(src).with_context(|| format!("failed to open {}", src.display()))?;
    let mut reader = BufReader::new(file);
    let mut writer = BufWriter::new(
        File::create(dst).with_context(|| format!("failed to create {}", dst.display()))?,
    );

    let mut names = HashSet::new();
    let mut complete = 0usize;
    let mut line = Vec::new();

    while complete < n {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }

        if let Some(rest) = line.strip_prefix(b"NAME") {
            names.insert(String::from_utf8_lossy(rest).trim().to_string());
        }

        writer.write_all(&line)?;

        if line.starts_with(b"//") {
            complete += 1;
        }
    }

    writer.flush()?;
    Ok(names)
}

/// Copy the stockholm records whose `#=GF ID` is in `names` from `src` to
/// `dst`, returning how many were found.
///
/// Stops once every name has been found, since Pfam is 500MB and the subsets
/// come off the front of it. On bytes for the same reason as [`subset_hmm`].
fn subset_sto(src: &Path, names: &HashSet<String>, dst: &Path) -> anyhow::Result<usize> {
    let file = File::open(src).with_context(|| format!("failed to open {}", src.display()))?;
    let mut reader = BufReader::new(file);
    let mut writer = BufWriter::new(
        File::create(dst).with_context(|| format!("failed to create {}", dst.display()))?,
    );

    let mut block: Vec<u8> = Vec::new();
    let mut id: Option<String> = None;
    let mut kept = 0usize;
    let mut line = Vec::new();

    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }

        if let Some(rest) = line.strip_prefix(b"#=GF ID") {
            id = Some(String::from_utf8_lossy(rest).trim().to_string());
        }

        block.extend_from_slice(&line);

        if line.starts_with(b"//") {
            if id.as_ref().is_some_and(|i| names.contains(i)) {
                writer.write_all(&block)?;
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
    Ok(kept)
}

/// Count the models in an hmm file by its `//` terminators, the same way the
/// subset does, so the ladder can be capped before any work happens.
fn count_models(path: &Path) -> anyhow::Result<usize> {
    let file =
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut reader = BufReader::new(file);

    let mut n = 0usize;
    let mut line = Vec::new();

    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.starts_with(b"//") {
            n += 1;
        }
    }

    Ok(n)
}

/// Draw every target rung in one pass, each a prefix of the next.
///
/// Returns the residues and bytes that landed in each.
fn deal(
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

struct Row {
    axis: &'static str,
    rung: usize,
    residues: u64,
    bytes: u64,
}

/// What each rung actually turned out to be, since sequences and models are not
/// uniform units of work and residues is the honest axis to plot against.
fn write_sizes(path: &Path, rows: &[Row]) -> anyhow::Result<()> {
    let headers = ["axis", "rung", "residues", "bytes"];

    let cells: Vec<[String; 4]> = rows
        .iter()
        .map(|r| {
            [
                r.axis.to_string(),
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
        for (c, &w) in row.iter().zip(&widths) {
            out.push_str(&pad(c, w));
            out.push(' ');
        }
        out.push('\n');
    }

    std::fs::write(path, out).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}
