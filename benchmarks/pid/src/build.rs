//! Assembles one benchmark: Pfam families split by identity, hidden in a
//! Swissprot decoy background.
//!
//! Everything external goes through the pipeline -- create-profmark, hmmbuild,
//! hmmemit -- so the build gets `--dry-run`, keeps the stderr of whatever
//! failed, and prints what it is doing while it does it. The assembly itself is
//! Rust, so it is a closure step in the same pipeline rather than something
//! that happens beside it.
//!
//! The profmark split is drawn once and shared: it depends only on Pfam and the
//! split parameters, and it is the expensive half.

use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, bail};
use clap::Parser;

use indexmap::IndexMap;
use libsail::collection::{Indexable, Iterable};
use libsail::seq::fasta::{DEFAULT_LINE_WIDTH, Fasta, FastaRecord};
use libsail::seq::stockholm::StockholmRecord;
use michi::{Closure, Cmd, PipelineBuilder, Progress, Step};

use rand::rngs::StdRng;
use rand::seq::{IndexedRandom, SliceRandom};
use rand::{RngExt, SeedableRng};

use crate::inputs;

/// Decoys per true pair in the target database.
const DECOY_RATIO: usize = 100;
/// Pairs above this percent identity are discarded; the benchmark targets the
/// twilight zone.
const PID_MAX: usize = 25;
/// Sequences shorter than this are unfindable at very low identity.
const SEQ_LEN_MIN: usize = 150;
const FAM_MIN: usize = 5;
const FAM_MAX: usize = 10;

/// The 20 standard residues, upper and lower case. Ambiguity codes (B, J, O,
/// U, X, Z) and gap characters deliberately do not count toward identity.
static AMINO: [bool; 256] = {
    let mut t = [false; 256];
    let letters = b"ACDEFGHIKLMNPQRSTVWY";
    let mut i = 0;
    while i < letters.len() {
        t[letters[i] as usize] = true;
        t[(letters[i] + 32) as usize] = true;
        i += 1;
    }
    t
};

#[derive(Parser, Debug)]
pub struct Args {
    /// Impose a limit on the number of benchmark pairs; the default keeps every
    /// pair that survives filtering. Decoys are added to the target database on
    /// top of this
    #[arg(short, long, value_name = "N")]
    pub pairs: Option<usize>,

    /// Seed for pair sampling and decoy generation, and for the profmark split
    #[arg(long, default_value_t = 67779, value_name = "N")]
    pub seed: u64,

    /// Maximum identity between the train and test halves of the split
    #[arg(long, default_value_t = 0.30, value_name = "X")]
    pub train_test_id: f64,

    /// Minimum test sequences per family
    #[arg(long, default_value_t = 10, value_name = "N")]
    pub min_test: usize,

    /// Maximum test sequences per family
    #[arg(long, default_value_t = 30, value_name = "N")]
    pub max_test: usize,

    /// Rebuild the profmark train/test split even if it already exists
    #[arg(long)]
    pub refresh_profmark: bool,

    /// Threads for hmmbuild
    #[arg(short, long, default_value_t = 8, value_name = "N")]
    pub threads: usize,

    #[arg(long)]
    pub dry_run: bool,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let src_sto = util::tools::pfam_sto()?;
    let src_fa = util::tools::swissprot()?;

    // resolved up front: the assembly takes long enough that finding out about
    // a missing hmmbuild afterwards would be miserable
    let hmmbuild = util::tools::hmmbuild()?;
    let hmmemit = util::tools::hmmemit()?;

    let pm = inputs::profmark();
    let split = args.refresh_profmark
        || !inputs::profmark_query().exists()
        || !inputs::profmark_target().exists();

    if !split {
        println!("reusing the profmark split in {}", pm.display());
    }

    let mut pl = PipelineBuilder::new().step(
        Cmd::new("mkdir")
            .name("dirs")
            .flag("-p")
            .path(&pm)
            .path(inputs::afa()),
    );

    if split {
        // --onlysplit names its output after the run, so the two halves are
        // renamed into the names everything downstream reads
        let stem = pm.join("benchmark");

        pl = pl
            .step(
                Step::serial([Cmd::new(util::tools::create_profmark()?)
                    .name("create-profmark")
                    .arg("-S", args.seed)
                    .arg("-1", format!("{:.2}", args.train_test_id))
                    .flag("--cluster")
                    .flag("--onlysplit")
                    .arg("--mintest", args.min_test)
                    .arg("--maxtest", args.max_test)
                    .path(&stem)
                    .path(&src_sto)])
                .name("profmark"),
            )
            .step(
                Step::from_closures([Closure::new("rename", move || {
                    fs::rename(stem.with_extension("train.msa"), inputs::profmark_query())?;
                    fs::rename(stem.with_extension("test.msa"), inputs::profmark_target())?;
                    fs::remove_file(stem.with_extension("tbl")).ok();
                    Ok(())
                })])
                .name("rename"),
            );
    }

    let pipeline = pl
        .step(
            Step::from_closures([Closure::new("assemble", {
                let (pairs, seed) = (args.pairs, args.seed);

                move || assemble(&src_sto, &src_fa, pairs, seed)
            })])
            .name("assemble"),
        )
        .step(
            Step::serial([Cmd::new(&hmmbuild)
                .name("hmmbuild")
                .arg("--cpu", args.threads)
                .path(inputs::query_hmm())
                .path(inputs::query_sto())])
            .name("profiles"),
        )
        .step(
            Step::serial([Cmd::new(&hmmemit)
                .name("hmmemit")
                .flag("-c")
                .path(inputs::query_hmm())
                .stdout_to(inputs::query_cons())])
            .name("consensus"),
        )
        .stderr_dir(inputs::tmp().join("build/stderr"))
        .sink(Progress::new())
        .build()
        .context("failed to build the assembly")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    // refused rather than cleared: there is one benchmark, assembling it is
    // expensive, and a rebuild in place would leave the last assembly's
    // families alongside this one's in afa/, which psiblast searches
    if inputs::exists() {
        bail!(
            "{} already exists; remove it to rebuild",
            inputs::dir().display()
        );
    }

    pipeline.run()?;

    println!("\nbuilt {}", inputs::dir().display());
    Ok(())
}

#[derive(Clone)]
struct Pair {
    pid: usize,
    family: String,
    query: String,
    target: String,
}

/// Assemble a benchmark from the profmark train/test split.
///
/// The RNG is seeded from the arguments rather than from entropy, so a given
/// pair limit and seed reproduce the same benchmark.
fn assemble(
    src_sto_path: &Path,
    src_fa_path: &Path,
    max_pairs: Option<usize>,
    seed: u64,
) -> anyhow::Result<()> {
    println!("loading alignments...");
    let mut query_sto =
        families(&inputs::profmark_query()).context("failed to parse the profmark query split")?;
    let mut target_sto = families(&inputs::profmark_target())
        .context("failed to parse the profmark target split")?;
    let src_sto = families(src_sto_path).context("failed to parse source sto")?;
    let src_fa = Fasta::open(src_fa_path).context("failed to parse source fasta")?;

    let afa_dir = inputs::afa();

    if target_sto.len() != query_sto.len() {
        bail!(
            "target/query family count mismatch: {} vs {}",
            target_sto.len(),
            query_sto.len()
        );
    }

    println!("{} sequence families found", target_sto.len());

    // short sequences at very low %ID are effectively impossible to find, so
    // they only add noise
    let len_filter = |fams: &mut IndexMap<String, StockholmRecord>| {
        for rec in fams.values_mut() {
            let keep = rec
                .names
                .iter()
                .zip(&rec.seqs)
                .filter(|(_, seq)| {
                    seq.iter().filter(|b| AMINO[**b as usize]).count() >= SEQ_LEN_MIN
                })
                .map(|(name, _)| name.clone())
                .collect::<HashSet<_>>();

            *rec = keep_rows(rec, &keep);
        }

        fams.retain(|_, rec| !rec.is_empty());
    };

    len_filter(&mut query_sto);
    len_filter(&mut target_sto);

    query_sto.retain(|_, rec| rec.depth() >= 10);

    // a family is only useful if it has both queries to search with and
    // targets to find
    target_sto.retain(|fam, _| query_sto.contains_key(fam));
    query_sto.retain(|fam, _| target_sto.contains_key(fam));

    println!(
        "{} families remain after length filter (>={SEQ_LEN_MIN})",
        target_sto.len()
    );

    // pair each target with its most similar query, so every target gets the
    // best shot it has
    let mut pairs_by_fam: HashMap<String, Vec<Pair>> =
        query_sto.keys().map(|fam| (fam.clone(), vec![])).collect();

    // one family's kept target rows, applied after the loop that reads them:
    // narrowing target_sto in place would need it borrowed both ways at once
    let mut kept: HashMap<String, HashSet<String>> = HashMap::new();

    for fam in query_sto.keys() {
        let src_rec = src_sto
            .get(fam)
            .with_context(|| format!("family {fam:?} missing from source stockholm"))?;

        let query_rec = query_sto.get(fam).expect("family present by construction");
        let target_rec = target_sto.get(fam).expect("family present by construction");

        let query_seqs = src_rec
            .names
            .iter()
            .zip(&src_rec.seqs)
            .filter(|(name, _)| query_rec.get(name).is_some())
            .collect::<Vec<_>>();

        let target_seqs = src_rec
            .names
            .iter()
            .zip(&src_rec.seqs)
            .filter(|(name, _)| target_rec.get(name).is_some())
            .collect::<Vec<_>>();

        let mut keep = HashSet::new();
        for (t_name, t_seq) in target_seqs.iter() {
            let mut best_pid = 0.0;
            let mut best_query = "";
            for (q_name, q_seq) in query_seqs.iter() {
                let pid = compute_pid(t_seq, q_seq);

                // profmark split the families by identity, so anything this
                // similar means the split did not do what we asked
                if pid > 0.5 {
                    bail!(
                        "unexpected {:.0}% identity between {t_name} and {q_name} in {fam}; \
                         check the profmark train/test split",
                        pid * 100.0
                    );
                }

                if pid > best_pid {
                    best_pid = pid;
                    best_query = q_name;
                }
            }

            let bin = (best_pid * 100.0).round() as usize;
            if bin <= PID_MAX {
                keep.insert((*t_name).clone());
                pairs_by_fam
                    .get_mut(fam)
                    .context("no pair vec for family")?
                    .push(Pair {
                        pid: bin,
                        family: fam.clone(),
                        query: best_query.to_string(),
                        target: t_name.to_string(),
                    })
            }
        }

        kept.insert(fam.clone(), keep);
    }

    for (fam, keep) in &kept {
        let rec = target_sto.get_mut(fam).context("family vanished")?;
        *rec = keep_rows(rec, keep);
    }

    pairs_by_fam.retain(|_, pairs| pairs.len() > FAM_MIN);
    pairs_by_fam.iter_mut().for_each(|(_, pairs)| {
        pairs.sort_by(|a, b| a.pid.cmp(&b.pid));
        pairs.truncate(FAM_MAX);
    });

    let mut rng = StdRng::seed_from_u64(seed);

    // families come out of a HashMap, so sort before sampling to keep the
    // seeded draw reproducible
    let mut fams: Vec<String> = pairs_by_fam.keys().cloned().collect();
    fams.sort();

    let mut pairs: Vec<Pair> = fams
        .iter()
        .flat_map(|f| pairs_by_fam.get(f).cloned().unwrap_or_default())
        .collect();

    pairs = match max_pairs {
        Some(max) if max < pairs.len() => pairs.sample(&mut rng, max).cloned().collect(),
        _ => pairs,
    };
    pairs.sort_by(|a, b| a.pid.cmp(&b.pid));

    println!("{} benchmark pairs", pairs.len());

    let mut tbl_writer = BufWriter::new(
        File::create(inputs::benchmark_tbl()).context("failed to open benchmark.tbl")?,
    );
    writeln!(tbl_writer, "#identity family target query")?;

    let mut targets: Vec<FastaRecord> = Vec::new();
    // a hash set because two targets can share a most-similar query
    let mut queries: HashSet<FastaRecord> = HashSet::new();

    let extract = |fams: &IndexMap<String, StockholmRecord>, fam: &str, seq: &str| {
        fams.get(fam)
            .and_then(|r| r.get(seq))
            .map(ungap)
            .context("failed to extract sequence from stockholm")
    };

    for (pair_idx, pair) in pairs.iter().enumerate() {
        let query = extract(&query_sto, &pair.family, &pair.query).map(|seq| FastaRecord {
            name: format!("{}|{}", pair.family, pair.query).into_bytes(),
            extra: Vec::new(),
            seq,
        })?;

        let target = extract(&target_sto, &pair.family, &pair.target).map(|seq| FastaRecord {
            name: format!("{}|{}|{}%:{}", pair.family, pair.target, pair.pid, pair_idx)
                .into_bytes(),
            extra: Vec::new(),
            seq,
        })?;

        targets.push(target);
        queries.insert(query);

        writeln!(
            tbl_writer,
            "{}% {} {} {}",
            pair.pid, pair.family, pair.target, pair.query
        )?;
    }

    let mut target_writer =
        BufWriter::new(File::create(inputs::target_fa()).context("failed to open target.fa")?);
    targets
        .iter()
        .try_for_each(|t| t.write_to(&mut target_writer, DEFAULT_LINE_WIDTH))?;

    let mut query_fa_writer =
        BufWriter::new(File::create(inputs::query_fa()).context("failed to open query.fa")?);

    let mut queries = queries.into_iter().collect::<Vec<_>>();
    queries.sort_by(|a, b| a.name.cmp(&b.name));
    queries
        .iter()
        .try_for_each(|q| q.write_to(&mut query_fa_writer, DEFAULT_LINE_WIDTH))?;

    let mut query_sto_writer =
        BufWriter::new(File::create(inputs::query_sto()).context("failed to open query.sto")?);
    fs::create_dir_all(&afa_dir)?;

    let query_names = queries
        .iter()
        .map(|q| -> anyhow::Result<&str> {
            Ok(q.name_str()?.split('|').next().expect("a split has a head"))
        })
        .collect::<anyhow::Result<HashSet<_>>>()?;

    query_sto.retain(|fam, _| query_names.contains(&fam.as_str()));

    query_sto
        .iter()
        .try_for_each(|(fam, rec)| -> anyhow::Result<()> {
            rec.write_to(&mut query_sto_writer)?;
            let mut afa_writer = BufWriter::new(File::create(afa_dir.join(format!("{fam}.afa")))?);

            let mut rows = rec.names.iter().zip(&rec.seqs);

            // the first sequence carries the family name so blast labels the
            // resulting profile usefully
            let (_, seq) = rows.next().expect("no seqs in query sto record");
            writeln!(afa_writer, ">{fam}")?;
            afa_writer.write_all(seq)?;
            writeln!(afa_writer)?;

            rows.try_for_each(|(name, seq)| -> anyhow::Result<()> {
                writeln!(afa_writer, ">{name}")?;
                afa_writer.write_all(seq)?;
                writeln!(afa_writer)?;
                Ok(())
            })
        })?;

    // ---- decoys ----

    let n_decoys = pairs.len() * DECOY_RATIO;
    println!(
        "sampling {n_decoys} decoys from {} source sequences...",
        src_fa.len()
    );

    // decoys are length-matched to real targets and shuffled, so they share the
    // benchmark's length and composition profile without any real homology
    let lengths: Vec<usize> = targets.iter().map(|t| t.seq.len()).collect();
    for decoy in decoys(&src_fa, &lengths, n_decoys, &mut rng)? {
        decoy.write_to(&mut target_writer, DEFAULT_LINE_WIDTH)?;
    }

    target_writer.flush()?;
    Ok(())
}

/// A Stockholm file keyed by family, which is what `#=GF ID` holds in
/// everything this benchmark reads.
fn families(path: &Path) -> anyhow::Result<IndexMap<String, StockholmRecord>> {
    // repaired in place rather than through from_utf8_lossy, which would hold
    // a second copy of the whole file
    let mut bytes =
        std::fs::read(path).with_context(|| format!("failed to read {}", path.display()))?;
    util::repair_utf8(&mut bytes);

    let sto = libsail::seq::stockholm::Stockholm::new(&bytes[..])
        .with_context(|| format!("failed to parse {}", path.display()))?;

    let mut out = IndexMap::with_capacity(sto.len());
    for rec in sto.iter() {
        let id = rec
            .id()
            .with_context(|| format!("a record in {} has no #=GF ID", path.display()))?;

        out.insert(id.to_string(), rec.clone());
    }

    Ok(out)
}

/// `rec` narrowed to the rows named in `keep`, carrying only their `#=GS`
/// lines.
///
/// A fresh record rather than an edit in place: `row` is the only thing that
/// keeps the name index true, and it only ever appends.
fn keep_rows(rec: &StockholmRecord, keep: &HashSet<String>) -> StockholmRecord {
    let mut out = StockholmRecord::default();

    out.gf = rec.gf.clone();
    out.gs = rec
        .gs
        .iter()
        .filter(|(name, _, _)| keep.contains(name))
        .cloned()
        .collect();

    for (name, seq) in rec.names.iter().zip(&rec.seqs) {
        if keep.contains(name) {
            let row = out.row(name);
            out.seqs[row] = seq.clone();
        }
    }

    out
}

/// An aligned row with its gap columns dropped.
fn ungap(row: &[u8]) -> Vec<u8> {
    row.iter()
        .copied()
        .filter(|&b| b != b'-' && b != b'.')
        .collect()
}

/// Generate decoy records by drawing a subsequence from `source` matched to the
/// length of a randomly chosen record in `lengths`, then shuffling it.
///
/// Shuffling preserves amino acid composition while destroying any real
/// homology, which is what makes a decoy a fair negative rather than simply an
/// unrelated sequence.
fn decoys(
    source: &Fasta,
    lengths: &[usize],
    count: usize,
    rng: &mut StdRng,
) -> anyhow::Result<Vec<FastaRecord>> {
    if source.is_empty() || lengths.is_empty() {
        bail!("cannot generate decoys from an empty source");
    }

    let mut out = Vec::with_capacity(count);
    let mut src: &[u8] = &[];

    for idx in 0..count {
        let decoy_len = lengths[rng.random_range(0..lengths.len())];

        // keep drawing until a source sequence is long enough to cut from
        while src.len() < decoy_len {
            src = &source
                .get(rng.random_range(0..source.len()))
                .context("bad source index")?
                .seq;
        }

        let start = rng.random_range(0..=src.len() - decoy_len);
        let mut seq = src[start..start + decoy_len].to_vec();
        seq.shuffle(rng);

        out.push(FastaRecord {
            name: format!("decoy{idx}").into_bytes(),
            extra: Vec::new(),
            seq,
        });

        src = &[];
    }

    Ok(out)
}

fn compute_pid(s1: &[u8], s2: &[u8]) -> f32 {
    debug_assert_eq!(s1.len(), s2.len());

    let mut match_cnt = 0usize;
    let mut pos_cnt = 0usize;

    s1.iter().zip(s2.iter()).for_each(|(&a, &b)| {
        if AMINO[a as usize] || AMINO[b as usize] {
            pos_cnt += 1;
            match_cnt += (a == b) as usize
        }
    });

    if pos_cnt == 0 {
        return 0.0;
    }

    match_cnt as f32 / pos_cnt as f32
}
