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

use bioio::fasta::{Fasta, FastaRecord};
use bioio::stockholm::Stockholm;
use pail::{Closure, Cmd, PipelineBuilder, Progress, Step};

use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::IndexedRandom;

use crate::inputs::{self, Inputs};

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
    /// Names the input set, under inputs/<size>/.
    #[arg(short, long, default_value = "toy")]
    pub size: String,

    /// Benchmark pairs to sample; 0 uses every pair that survives filtering.
    /// Decoys are added to the target database on top of this.
    #[arg(short, long, default_value_t = 50)]
    pub pairs: usize,

    /// Seed for pair sampling and decoy generation, and for the profmark split.
    #[arg(long, default_value_t = 67779)]
    pub seed: u64,

    /// Maximum identity between the train and test halves of the split.
    #[arg(long, default_value_t = 0.30)]
    pub train_test_id: f64,

    /// Minimum test sequences per family.
    #[arg(long, default_value_t = 10)]
    pub min_test: usize,

    /// Maximum test sequences per family.
    #[arg(long, default_value_t = 30)]
    pub max_test: usize,

    /// Rebuild the profmark train/test split even if it already exists.
    #[arg(long)]
    pub refresh_profmark: bool,

    /// Threads for hmmbuild.
    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,

    #[arg(long)]
    pub dry_run: bool,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let src_sto = tools::pfam_sto()?;
    let src_fa = tools::swissprot()?;

    // resolved up front: the assembly takes long enough that finding out about
    // a missing hmmbuild afterwards would be miserable
    let hmmbuild = tools::hmmbuild()?;
    let hmmemit = tools::hmmemit()?;

    let set = Inputs::new(&args.size);

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
            .path(set.afa()),
    );

    if split {
        // --onlysplit names its output after the run, so the two halves are
        // renamed into the names everything downstream reads
        let stem = pm.join("benchmark");

        pl = pl
            .step(
                Step::serial([Cmd::new(tools::create_profmark()?)
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
                let set = Inputs::new(&args.size);
                let (pairs, seed) = (args.pairs, args.seed);

                move || assemble(&set, &src_sto, &src_fa, (pairs > 0).then_some(pairs), seed)
            })])
            .name("assemble"),
        )
        .step(
            Step::serial([Cmd::new(&hmmbuild)
                .name("hmmbuild")
                .arg("--cpu", args.threads)
                .path(set.query_hmm())
                .path(set.query_sto())])
            .name("profiles"),
        )
        .step(
            Step::serial([Cmd::new(&hmmemit)
                .name("hmmemit")
                .flag("-c")
                .path(set.query_hmm())
                .stdout_to(set.query_cons())])
            .name("consensus"),
        )
        .stderr_dir(set.run_dir().join("tmp/stderr"))
        .sink(Progress::new())
        .build()
        .context("failed to build the assembly")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    // a rebuild in place would leave the last assembly's families alongside
    // this one's in afa/, and psiblast searches the directory
    if set.exists() {
        fs::remove_dir_all(set.dir())
            .with_context(|| format!("failed to clear {}", set.dir().display()))?;
    }

    pipeline.run()?;

    println!("\nbuilt {}", set.dir().display());
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
/// size and seed reproduce the same benchmark.
fn assemble(
    set: &Inputs,
    src_sto_path: &Path,
    src_fa_path: &Path,
    max_pairs: Option<usize>,
    seed: u64,
) -> anyhow::Result<()> {
    println!("loading alignments...");
    let mut query_sto = Stockholm::from_path(inputs::profmark_query())
        .context("failed to parse the profmark query split")?;
    let mut target_sto = Stockholm::from_path(inputs::profmark_target())
        .context("failed to parse the profmark target split")?;
    let src_sto = Stockholm::from_path(src_sto_path).context("failed to parse source sto")?;
    let src_fa = Fasta::from_path(src_fa_path).context("failed to parse source fasta")?;

    let afa_dir = set.afa();

    if target_sto.records.len() != query_sto.records.len() {
        bail!(
            "target/query family count mismatch: {} vs {}",
            target_sto.records.len(),
            query_sto.records.len()
        );
    }

    println!("{} sequence families found", target_sto.records.len());

    // short sequences at very low %ID are effectively impossible to find, so
    // they only add noise
    let len_filter = |sto: &mut Stockholm| {
        sto.records.values_mut().for_each(|rec| {
            let keep = rec
                .sequences
                .iter()
                .filter(|(_, seq)| {
                    seq.bytes().filter(|b| AMINO[*b as usize]).count() >= SEQ_LEN_MIN
                })
                .map(|(name, _)| name.clone())
                .collect::<HashSet<_>>();

            rec.gs_meta.retain(|name, _| keep.contains(name));
            rec.sequences.retain(|name, _| keep.contains(name));
        });

        sto.records.retain(|_, rec| !rec.sequences.is_empty());
    };

    len_filter(&mut query_sto);
    len_filter(&mut target_sto);

    query_sto.records.retain(|_, rec| rec.sequences.len() >= 10);

    // a family is only useful if it has both queries to search with and
    // targets to find
    target_sto
        .records
        .retain(|fam, _| query_sto.records.contains_key(fam));
    query_sto
        .records
        .retain(|fam, _| target_sto.records.contains_key(fam));

    println!(
        "{} families remain after length filter (>={SEQ_LEN_MIN})",
        target_sto.records.len()
    );

    // pair each target with its most similar query, so every target gets the
    // best shot it has
    let mut pairs_by_fam: HashMap<String, Vec<Pair>> = query_sto
        .records
        .keys()
        .map(|fam| (fam.clone(), vec![]))
        .collect();

    for fam in query_sto.records.keys() {
        let src_seqs = &src_sto
            .get(fam)
            .map(|rec| &rec.sequences)
            .with_context(|| format!("family {fam:?} missing from source stockholm"))?;

        let query_seqs = src_seqs
            .iter()
            .filter(|(name, _)| {
                query_sto
                    .get(fam)
                    .expect("family present by construction")
                    .sequences
                    .contains_key(*name)
            })
            .collect::<Vec<_>>();

        let target_seqs = src_seqs
            .iter()
            .filter(|(name, _)| {
                target_sto
                    .get(fam)
                    .expect("family present by construction")
                    .sequences
                    .contains_key(*name)
            })
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
                keep.insert(t_name);
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

        let sto = target_sto.get_mut(fam).context("family vanished")?;
        sto.sequences.retain(|n, _| keep.contains(&n));
        sto.gs_meta.retain(|n, _| keep.contains(&n));
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
        Some(max) if max < pairs.len() => pairs.choose_multiple(&mut rng, max).cloned().collect(),
        _ => pairs,
    };
    pairs.sort_by(|a, b| a.pid.cmp(&b.pid));

    println!("{} benchmark pairs", pairs.len());

    let mut tbl_writer =
        BufWriter::new(File::create(set.benchmark_tbl()).context("failed to open benchmark.tbl")?);
    writeln!(tbl_writer, "#identity family target query")?;

    let mut targets: Vec<FastaRecord> = Vec::new();
    // a hash set because two targets can share a most-similar query
    let mut queries: HashSet<FastaRecord> = HashSet::new();

    let extract = |sto: &Stockholm, fam: &str, seq: &str| {
        sto.get(fam)
            .and_then(|r| r.get(seq))
            .map(|s| s.replace(['-', '.'], ""))
            .context("failed to extract sequence from stockholm")
    };

    for (pair_idx, pair) in pairs.iter().enumerate() {
        let query = extract(&query_sto, &pair.family, &pair.query).map(|seq| FastaRecord {
            name: format!("{}|{}", pair.family, pair.query),
            extra: String::new(),
            seq,
        })?;

        let target = extract(&target_sto, &pair.family, &pair.target).map(|seq| FastaRecord {
            name: format!("{}|{}|{}%:{}", pair.family, pair.target, pair.pid, pair_idx),
            extra: String::new(),
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
        BufWriter::new(File::create(set.target_fa()).context("failed to open target.fa")?);
    targets
        .iter()
        .try_for_each(|t| writeln!(target_writer, "{t}"))?;

    let mut query_fa_writer =
        BufWriter::new(File::create(set.query_fa()).context("failed to open query.fa")?);

    let mut queries = queries.into_iter().collect::<Vec<_>>();
    queries.sort_by(|a, b| a.name.cmp(&b.name));
    queries
        .iter()
        .try_for_each(|q| writeln!(query_fa_writer, "{q}"))?;

    let mut query_sto_writer =
        BufWriter::new(File::create(set.query_sto()).context("failed to open query.sto")?);
    fs::create_dir_all(&afa_dir)?;

    let query_names = queries
        .iter()
        .map(|q| q.name.split('|').next().unwrap())
        .collect::<HashSet<_>>();

    query_sto
        .records
        .retain(|fam, _| query_names.contains(&fam.as_str()));

    query_sto.records.iter().try_for_each(|(fam, rec)| {
        writeln!(query_sto_writer, "{rec}")?;
        let mut afa_writer = BufWriter::new(File::create(afa_dir.join(format!("{fam}.afa")))?);

        let mut seqs = rec.sequences.iter();

        // the first sequence carries the family name so blast labels the
        // resulting profile usefully
        let (_, seq) = seqs.next().expect("no seqs in query sto record");
        writeln!(afa_writer, ">{fam}")?;
        writeln!(afa_writer, "{seq}")?;

        seqs.try_for_each(|(name, seq)| {
            writeln!(afa_writer, ">{name}")?;
            writeln!(afa_writer, "{seq}")
        })
    })?;

    // ---- decoys ----

    let n_decoys = pairs.len() * DECOY_RATIO;
    println!(
        "sampling {n_decoys} decoys from {} source sequences...",
        src_fa.records.len()
    );

    // decoys are length-matched to real targets and shuffled, so they share the
    // benchmark's length and composition profile without any real homology
    let lengths: Vec<usize> = targets.iter().map(|t| t.seq.len()).collect();
    for decoy in bioio::fasta::decoys(&src_fa, &lengths, n_decoys, &mut rng)? {
        writeln!(target_writer, "{decoy}")?;
    }

    target_writer.flush()?;
    Ok(())
}

fn compute_pid(s1: &str, s2: &str) -> f32 {
    debug_assert_eq!(s1.len(), s2.len());

    let mut match_cnt = 0usize;
    let mut pos_cnt = 0usize;

    s1.as_bytes()
        .iter()
        .zip(s2.as_bytes().iter())
        .for_each(|(&a, &b)| {
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
