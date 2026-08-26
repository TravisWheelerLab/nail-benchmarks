//! Learning per-family false-positive score cutoffs from reversed decoys.
//!
//! Reversing a protein keeps its composition but destroys its homology, so a
//! reversed sequence that still scores against a family is measuring noise.
//! Collecting those scores per family gives a threshold above which a hit is
//! unlikely to be chance.
//!
//! Doing that exhaustively would mean searching every family against every
//! reversed sequence. Instead it runs in two stages:
//!
//!   1. `recruit` — a cheap sweep of every family against the reversed shards,
//!      which finds the small subset of sequences that score at all.
//!   2. `search` — an exhaustive pass, one family at a time, against just its
//!      own recruits, in both directions.
//!
//! Stage 2 is the one the cutoffs are read from. It runs at high sensitivity
//! with the prefilter effectively disabled, so a decoy's score is its real
//! score rather than one truncated by stage 1's parameters. It also searches
//! the *forward* sequences, so a recruit that turns out to be a genuine family
//! member can be dropped instead of inflating the threshold.
//!
//! `decoys` sits between them: it reads stage 1's hit tables, un-reverses the
//! sequences that hit, and splits the query set per family.
//!
//! `recruit` is one big search per shard, so it runs through `pipeline` the
//! same way the rest of this crate does. `search` is the opposite shape — many
//! tiny per-family searches, run several at once rather than one at a time —
//! which `pipeline`'s `Cmd`/`Step` DAG has no batched form for (a family needs
//! several commands in a fixed order: build its query profile, then search
//! and convert per direction). It runs as a plain rayon pool instead, the same
//! way `reverse`, `decoys` and `learn` already do their own parallel work in
//! this file.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{Context, bail};
use clap::{Parser, Subcommand};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use bioio::tbl::{BlastTable, HitTable, HmmerTable, NailTable};
use pail::{Cmd as PCmd, PipelineBuilder, Progress, Step};
use tools::{hmmsearch, mmseqs, nail};

use crate::inputs::{self, shards};

/// Name of the calibration directory when none is given.
pub const DEFAULT_OUT: &str = "cutoffs";

// every stage reports down to here, so scores are comparable
const EVALUE: &str = "10";

// recruitment only has to nominate candidates, so it runs a cheap sweep
// rather than the decoy stage's wide-open one
const RECRUIT_S: &str = "11.0";
const RECRUIT_MAX_SEQS: &str = "5000";

// wide open, so a decoy's score is its real score rather than one truncated
// by a prefilter. one thread each: the parallelism is in running many
// families at once, not in any one search
const DECOY_S: &str = "12.0";
const DECOY_MAX_SEQS: &str = "1000000000";

const NAIL_RECRUIT: &str = "nail-recruit";
const MMSEQS_RECRUIT: &str = "mmseqs-recruit";

const NAIL_DECOY: &str = "nail";
const MMSEQS_DECOY: &str = "mmseqs";
const HMMER_DECOY: &str = "hmmer";

// ------------------------------------------------------------------ layout

/// Where every artifact of a calibration lives.
///
/// The whole thing hangs off one directory at the crate root, so several
/// calibrations — different shard counts, different parameters — can sit side
/// by side and be deleted as a unit. The inputs it reads are the shared ones
/// every pipeline reads.
struct Layout {
    root: PathBuf,
}

impl Layout {
    fn new(out: &str) -> anyhow::Result<Self> {
        for dir in [inputs::fixed::queries(), inputs::fixed::targets()] {
            if !dir.is_dir() {
                bail!("{} does not exist; run `mgy build` first", dir.display());
            }
        }

        Ok(Layout {
            root: crate::dir().join(out),
        })
    }

    /// Forward shards, built by `mgy build fixed`.
    fn targets(&self) -> PathBuf {
        inputs::fixed::targets()
    }

    fn query_hmm(&self) -> PathBuf {
        inputs::fixed::query_hmm()
    }

    fn query_sto(&self) -> PathBuf {
        inputs::fixed::query_sto()
    }

    /// The mmseqs profile db `mgy build fixed` made from the query set. Reused here
    /// rather than rebuilt, since recruit searches the whole query set.
    fn query_db(&self) -> PathBuf {
        inputs::fixed::query_db()
    }

    fn targets_rev(&self) -> PathBuf {
        self.root.join("targets-rev")
    }

    fn recruit(&self) -> PathBuf {
        self.root.join("recruit")
    }

    fn recruit_results(&self) -> PathBuf {
        self.recruit().join("results")
    }

    fn decoys(&self) -> PathBuf {
        self.root.join("decoys")
    }

    fn decoys_rev(&self) -> PathBuf {
        self.root.join("decoys-rev")
    }

    fn queries_hmm(&self) -> PathBuf {
        self.root.join("queries/hmm")
    }

    fn queries_sto(&self) -> PathBuf {
        self.root.join("queries/sto")
    }

    fn results(&self) -> PathBuf {
        self.root.join("results")
    }

    fn cutoffs_txt(&self) -> PathBuf {
        self.root.join("cutoffs.txt")
    }
}

// ---------------------------------------------------------------- commands

#[derive(Subcommand)]
pub enum Cmd {
    /// Reverse the target shards into the calibration directory.
    Reverse(ReverseArgs),
    /// Sweep every family against the reversed shards to find candidates.
    Recruit(RecruitArgs),
    /// Un-reverse what hit, group it per family, and split the query set.
    Decoys(DecoysArgs),
    /// Search each family against its own decoys, forward and reversed.
    Search(SearchArgs),
    /// Turn the decoy scores into per-family cutoffs.
    Learn(LearnArgs),
    /// Run every stage in order.
    All(AllArgs),
}

/// Which calibration directory to work in.
#[derive(Parser, Debug, Clone)]
pub struct Where {
    /// Calibration directory, created under benchmarks/mgy/.
    #[arg(long, default_value = DEFAULT_OUT)]
    pub out: String,
}

#[derive(Parser, Debug)]
pub struct ReverseArgs {
    #[command(flatten)]
    pub place: Where,

    /// Reverse only the first N shards. This is the one place the size of the
    /// calibration set is decided; every later stage uses whatever is here.
    #[arg(short = 'n', long)]
    pub shards: Option<usize>,

    /// Threads for the reversal itself.
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Parser, Debug)]
pub struct RecruitArgs {
    #[command(flatten)]
    pub place: Where,

    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,

    /// List the commands and exit without executing anything.
    #[arg(long)]
    pub dry_run: bool,
}

#[derive(Parser, Debug)]
pub struct DecoysArgs {
    #[command(flatten)]
    pub place: Where,

    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Parser, Debug)]
pub struct SearchArgs {
    #[command(flatten)]
    pub place: Where,

    /// How many families to search at once. Each search is single-threaded,
    /// so this is the whole of the parallelism.
    #[arg(short = 'j', long)]
    pub jobs: Option<usize>,
}

#[derive(Parser, Debug)]
pub struct LearnArgs {
    #[command(flatten)]
    pub place: Where,

    /// Forward hits at or below this E-value are treated as real, and their
    /// reversed counterparts are excluded from the decoy scores.
    #[arg(short = 'e', default_value_t = 1e-3, value_name = "F")]
    pub reverse_e_cutoff: f64,

    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Parser, Debug)]
pub struct AllArgs {
    #[command(flatten)]
    pub place: Where,

    #[arg(short = 'n', long)]
    pub shards: Option<usize>,

    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,

    #[arg(short = 'j', long)]
    pub jobs: Option<usize>,
}

pub fn main(cmd: Cmd) -> anyhow::Result<()> {
    match cmd {
        Cmd::Reverse(args) => reverse(args),
        Cmd::Recruit(args) => recruit(args),
        Cmd::Decoys(args) => decoys(args),
        Cmd::Search(args) => search(args),
        Cmd::Learn(args) => learn(args),
        Cmd::All(args) => all(args),
    }
}

// ----------------------------------------------------------------- reverse

fn reverse(args: ReverseArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.out)?;
    let src = layout.targets();
    let dst = layout.targets_rev();

    let mut found = shards(&src)?;
    if let Some(n) = args.shards {
        if n == 0 {
            bail!("--shards 0 would leave nothing to calibrate against");
        }
        if n > found.len() {
            eprintln!(
                "warning: asked for {n} shards but {} has only {}; using all of them",
                src.display(),
                found.len()
            );
        }
        found.truncate(n);
    }

    if dst.exists() {
        std::fs::remove_dir_all(&dst)?;
    }
    std::fs::create_dir_all(&dst)?;

    println!("reversing {} shards into {}...", found.len(), dst.display());

    let pool = pool(args.threads)?;
    pool.install(|| {
        found
            .par_iter()
            .try_for_each(|(i, path)| -> anyhow::Result<()> {
                // plain `<n>.fa`, not `<n>.rev.fa`: the directory already says
                // these are reversed, and the shard index has to stay readable off
                // the stem for the stages downstream
                bioio::fasta::reverse(path, &dst.join(format!("{i}.fa")))
                    .with_context(|| format!("failed to reverse shard {i}"))
            })
    })?;

    println!("reversed {} shards", found.len());
    Ok(())
}

// ----------------------------------------------------------------- recruit

fn recruit(args: RecruitArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.out)?;
    let rev = layout.targets_rev();

    if !rev.is_dir() {
        bail!(
            "no reversed shards in {}; run `mgy cutoffs reverse` first",
            rev.display()
        );
    }

    let nail_bin = nail()?;
    let mmseqs_bin = mmseqs()?;

    let query_hmm = layout.query_hmm();
    let query_db = layout.query_db();

    let recruit_dir = layout.recruit();
    let results = layout.recruit_results();
    let tmp = recruit_dir.join("tmp");

    let mut pl = PipelineBuilder::new().step(PCmd::new("mkdir").flag("-p").path(&results));

    for (idx, shard) in shards(&rev)? {
        let scratch = tmp.join(format!("shard-{idx}"));
        let target_db = scratch.join("targetDB/targetDB");
        let aln_db = scratch.join("alnDB/alnDB");

        pl = pl
            .step(
                Step::serial([
                    PCmd::new("mkdir")
                        .name("dirs")
                        .flag("-p")
                        .path(scratch.join("targetDB"))
                        .path(scratch.join("alnDB")),
                    PCmd::new(&mmseqs_bin)
                        .name("createdb")
                        .sub("createdb")
                        .path(&shard)
                        .path(&target_db),
                ])
                .name(format!("prep.{idx}")),
            )
            .step(
                Step::serial([PCmd::new(&nail_bin)
                    .sub("search")
                    .arg("--mmseqs-path", &mmseqs_bin)
                    .arg("-t", args.threads)
                    .arg("--tmp-dir", scratch.join("nail"))
                    .arg("--mmseqs-s", RECRUIT_S)
                    .arg("--mmseqs-max-seqs", RECRUIT_MAX_SEQS)
                    .arg("-E", EVALUE)
                    .arg(
                        "--tbl-out",
                        results.join(format!("{NAIL_RECRUIT}.{idx}.tbl")),
                    )
                    .flag("--allow-overwrite")
                    .path(&query_hmm)
                    .path(&shard)])
                .name(format!("nail.{idx}")),
            )
            .step(
                Step::serial([
                    PCmd::new(&mmseqs_bin)
                        .name("search")
                        .sub("search")
                        .arg("--threads", args.threads)
                        .arg("-s", RECRUIT_S)
                        .arg("--max-seqs", RECRUIT_MAX_SEQS)
                        .arg("-e", EVALUE)
                        .path(&query_db)
                        .path(&target_db)
                        .path(&aln_db)
                        .path(scratch.join("work")),
                    PCmd::new(&mmseqs_bin)
                        .name("convertalis")
                        .sub("convertalis")
                        .arg("--format-mode", 0)
                        .path(&query_db)
                        .path(&target_db)
                        .path(&aln_db)
                        .path(results.join(format!("{MMSEQS_RECRUIT}.{idx}.tbl"))),
                ])
                .name(format!("mmseqs.{idx}")),
            )
            .step(PCmd::new("rm").name("clean").flag("-rf").path(&scratch));
    }

    let pipeline = pl
        .stderr_dir(tmp.join("stderr"))
        .sink(Progress::new())
        .build()?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    pipeline.run()
}

// ------------------------------------------------------------------ decoys

fn decoys(args: DecoysArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.out)?;
    let recruit_results = layout.recruit_results();

    let shard_list: Vec<String> = shards(&layout.targets_rev())?
        .into_iter()
        .map(|(i, _)| i.to_string())
        .collect();

    println!("reading {} recruited shards...", shard_list.len());

    let pool = pool(args.threads)?;

    // per shard, which families each recruited sequence hit. Only names travel
    // here; the sequences themselves are read in the second pass.
    let wanted: Vec<(String, HashMap<String, Vec<String>>)> = pool.install(|| {
        shard_list
            .par_iter()
            .map(
                |shard| -> anyhow::Result<(String, HashMap<String, Vec<String>>)> {
                    let mut map: HashMap<String, Vec<String>> = HashMap::new();

                    let path = recruit_results.join(format!("{NAIL_RECRUIT}.{shard}.tbl"));
                    let tbl = HitTable::from_path::<_, NailTable>(&path)
                        .with_context(|| format!("failed to read {}", path.display()))?;
                    collect(tbl, &mut map);

                    let path = recruit_results.join(format!("{MMSEQS_RECRUIT}.{shard}.tbl"));
                    let tbl = HitTable::from_path::<_, BlastTable>(&path)
                        .with_context(|| format!("failed to read {}", path.display()))?;
                    collect(tbl, &mut map);

                    for v in map.values_mut() {
                        v.sort();
                        v.dedup();
                    }

                    Ok((shard.clone(), map))
                },
            )
            .collect::<anyhow::Result<Vec<_>>>()
    })?;

    let families: HashSet<String> = wanted.iter().flat_map(|(_, m)| m.keys().cloned()).collect();

    if families.is_empty() {
        bail!("no family recruited any decoys; is the recruit stage's E-value too strict?");
    }

    println!("{} families recruited decoys", families.len());

    // ---- un-reverse: pull the forward sequences the reversed ones came from

    let decoy_dir = layout.decoys();
    for dir in [&decoy_dir, &layout.decoys_rev()] {
        if dir.exists() {
            std::fs::remove_dir_all(dir)?;
        }
    }
    std::fs::create_dir_all(&decoy_dir)?;

    // one lock per family: shards are read in parallel and any of them may
    // contribute to any family
    let handles: HashMap<&str, Mutex<PathBuf>> = families
        .iter()
        .map(|f| (f.as_str(), Mutex::new(decoy_dir.join(format!("{f}.fa")))))
        .collect();

    let targets = layout.targets();
    pool.install(|| {
        wanted
            .par_iter()
            .try_for_each(|(shard, map)| -> anyhow::Result<()> {
                // invert to a name lookup so the shard is read once, streaming,
                // rather than held in memory
                let mut by_target: HashMap<&str, Vec<&str>> = HashMap::new();
                for (family, names) in map {
                    for name in names {
                        by_target
                            .entry(name.as_str())
                            .or_default()
                            .push(family.as_str());
                    }
                }

                let path = targets.join(format!("{shard}.fa"));
                let mut reader = bioio::fasta::Reader::from_path(&path)
                    .with_context(|| format!("failed to open {}", path.display()))?;

                let mut buffers: HashMap<&str, String> = HashMap::new();
                while let Some(rec) = reader.next_record()? {
                    let Some(fams) = by_target.get(rec.name.as_str()) else {
                        continue;
                    };
                    for family in fams {
                        let buf = buffers.entry(family).or_default();
                        buf.push_str(&format!("{rec}\n"));
                    }
                }

                for (family, text) in buffers {
                    let guard = handles
                        .get(family)
                        .with_context(|| format!("no handle for family {family}"))?
                        .lock()
                        .expect("family mutex poisoned");

                    let mut file = std::fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(&*guard)?;
                    file.write_all(text.as_bytes())?;
                }

                Ok(())
            })
    })?;

    // ---- reverse them back, so both directions are searched in stage 2

    let rev_dir = layout.decoys_rev();
    std::fs::create_dir_all(&rev_dir)?;

    println!("reversing decoys...");
    let names: Vec<&String> = families.iter().collect();
    pool.install(|| {
        names.par_iter().try_for_each(|f| -> anyhow::Result<()> {
            bioio::fasta::reverse(
                &decoy_dir.join(format!("{f}.fa")),
                &rev_dir.join(format!("{f}.rev.fa")),
            )
            .with_context(|| format!("failed to reverse decoys for {f}"))
        })
    })?;

    // ---- split the query set, so each family can be searched on its own

    println!("splitting queries for {} families...", families.len());

    let hmm = bioio::hmm::explode(layout.query_hmm(), &families, layout.queries_hmm())?;
    let sto = bioio::stockholm::explode(layout.query_sto(), &families, layout.queries_sto())?;

    if hmm != families.len() || sto != families.len() {
        bail!(
            "recruited {} families but found {hmm} in query.hmm and {sto} in query.sto",
            families.len()
        );
    }

    println!("wrote {}", layout.decoys().display());
    Ok(())
}

/// Fold a hit table into a family to target-name map.
fn collect(tbl: HitTable, map: &mut HashMap<String, Vec<String>>) {
    for (query, hits) in tbl.to_query_map() {
        let entry = map.entry(query).or_default();
        entry.extend(hits.iter().map(|h| h.target.clone()));
    }
}

// ------------------------------------------------------------------ search

/// Run one command to completion, its combined stdout and stderr appended to
/// `log`. mmseqs reports failures on stdout rather than stderr, so both share
/// one file rather than splitting a failure across two.
fn run(cmd: &mut Command, log: &Path) -> anyhow::Result<()> {
    let out = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(log)
        .with_context(|| format!("failed to open {}", log.display()))?;

    let status = cmd
        .stdout(out.try_clone()?)
        .stderr(out)
        .status()
        .with_context(|| format!("failed to spawn {cmd:?}"))?;

    if !status.success() {
        bail!("{cmd:?} exited with {status}; see {}", log.display());
    }

    Ok(())
}

fn search(args: SearchArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.out)?;
    let decoy_dir = layout.decoys();

    if !decoy_dir.is_dir() {
        bail!(
            "no decoys in {}; run `mgy cutoffs decoys` first",
            decoy_dir.display()
        );
    }

    let mut families: Vec<String> = std::fs::read_dir(&decoy_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "fa"))
        .filter_map(|p| Some(p.file_stem()?.to_str()?.to_string()))
        .collect();
    families.sort();

    if families.is_empty() {
        bail!("no decoy files in {}", decoy_dir.display());
    }

    let rev_dir = layout.decoys_rev();
    let hmm_dir = layout.queries_hmm();
    let sto_dir = layout.queries_sto();

    let results = layout.results();
    if results.exists() {
        std::fs::remove_dir_all(&results)?;
    }
    std::fs::create_dir_all(&results)?;

    let tmp = layout.root.join("tmp");
    std::fs::create_dir_all(&tmp)?;

    let nail_bin = nail()?;
    let mmseqs_bin = mmseqs()?;
    let hmmsearch_bin = hmmsearch()?;

    let jobs = args
        .jobs
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(4, |n| n.get()));

    println!("searching {} families, {jobs} at once...", families.len());

    let done = AtomicUsize::new(0);
    let total = families.len();

    let pool = pool(jobs)?;
    pool.install(|| {
        families
            .par_iter()
            .try_for_each(|family| -> anyhow::Result<()> {
                let scratch = tmp.join(family);
                std::fs::create_dir_all(&scratch)?;
                let log = scratch.join("log");

                let hmm = hmm_dir.join(format!("{family}.hmm"));
                let sto = sto_dir.join(format!("{family}.sto"));

                // the query profile is the same for both directions, so it is
                // built once per family rather than once per search
                let msa_db = scratch.join("msaDB");
                let query_db = scratch.join("queryDB");
                run(
                    Command::new(&mmseqs_bin)
                        .arg("convertmsa")
                        .arg(&sto)
                        .arg(&msa_db)
                        .arg("--identifier-field")
                        .arg("0"),
                    &log,
                )?;
                run(
                    Command::new(&mmseqs_bin)
                        .arg("msa2profile")
                        .arg(&msa_db)
                        .arg(&query_db)
                        .arg("--match-mode")
                        .arg("1"),
                    &log,
                )?;

                let directions = [
                    (family.clone(), decoy_dir.join(format!("{family}.fa"))),
                    (
                        format!("{family}.rev"),
                        rev_dir.join(format!("{family}.rev.fa")),
                    ),
                ];

                for (label, target) in directions {
                    let dir_scratch = scratch.join(&label);
                    std::fs::create_dir_all(&dir_scratch)?;

                    run(
                        Command::new(&nail_bin)
                            .arg("search")
                            .arg("--mmseqs-path")
                            .arg(&mmseqs_bin)
                            .arg("-t")
                            .arg("1")
                            .arg("--tmp-dir")
                            .arg(dir_scratch.join("nail"))
                            .arg("--mmseqs-s")
                            .arg(DECOY_S)
                            .arg("--mmseqs-max-seqs")
                            .arg(DECOY_MAX_SEQS)
                            .arg("-E")
                            .arg(EVALUE)
                            .arg("--allow-overwrite")
                            .arg("--tbl-out")
                            .arg(results.join(format!("{NAIL_DECOY}.{label}.tbl")))
                            .arg(&hmm)
                            .arg(&target),
                        &log,
                    )?;

                    let target_db = dir_scratch.join("targetDB/targetDB");
                    let aln_db = dir_scratch.join("alnDB/alnDB");
                    std::fs::create_dir_all(target_db.parent().unwrap())?;
                    std::fs::create_dir_all(aln_db.parent().unwrap())?;

                    run(
                        Command::new(&mmseqs_bin)
                            .arg("createdb")
                            .arg(&target)
                            .arg(&target_db),
                        &log,
                    )?;
                    run(
                        Command::new(&mmseqs_bin)
                            .arg("search")
                            .arg(&query_db)
                            .arg(&target_db)
                            .arg(&aln_db)
                            .arg(dir_scratch.join("work"))
                            .arg("--threads")
                            .arg("1")
                            .arg("-s")
                            .arg(DECOY_S)
                            .arg("--max-seqs")
                            .arg(DECOY_MAX_SEQS)
                            .arg("-e")
                            .arg(EVALUE),
                        &log,
                    )?;
                    run(
                        Command::new(&mmseqs_bin)
                            .arg("convertalis")
                            .arg(&query_db)
                            .arg(&target_db)
                            .arg(&aln_db)
                            .arg(results.join(format!("{MMSEQS_DECOY}.{label}.tbl")))
                            .arg("--format-mode")
                            .arg("0"),
                        &log,
                    )?;

                    run(
                        Command::new(&hmmsearch_bin)
                            .arg("--cpu")
                            .arg("1")
                            .arg("-E")
                            .arg(EVALUE)
                            .arg("-o")
                            .arg("/dev/null")
                            .arg("--tblout")
                            .arg(results.join(format!("{HMMER_DECOY}.{label}.tbl")))
                            .arg("--domtblout")
                            .arg(results.join(format!("{HMMER_DECOY}.{label}.domtbl")))
                            .arg(&hmm)
                            .arg(&target),
                        &log,
                    )?;
                }

                std::fs::remove_dir_all(&scratch).ok();

                let n = done.fetch_add(1, Ordering::Relaxed) + 1;
                if n.is_multiple_of(50) || n == total {
                    eprint!("\r  {n}/{total} families searched");
                }

                Ok(())
            })
    })?;
    eprintln!();

    println!("wrote {}", results.display());
    Ok(())
}

// ------------------------------------------------------------------- learn

/// How many decoy scores are recorded per family per tool.
const N_SCORES: usize = 5;

fn learn(args: LearnArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.out)?;
    let results = layout.results();

    let mut families: Vec<String> = std::fs::read_dir(layout.decoys())?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "fa"))
        .filter_map(|p| Some(p.file_stem()?.to_str()?.to_string()))
        .collect();
    families.sort();

    let out_path = layout.cutoffs_txt();
    let out = Mutex::new(BufWriter::new(File::create(&out_path)?));
    let skipped = AtomicUsize::new(0);

    let pool = pool(args.threads)?;
    pool.install(|| {
        families
            .par_iter()
            .try_for_each(|family| -> anyhow::Result<()> {
                let nail_scores =
                    decoy_scores::<NailTable>(&results, NAIL_DECOY, family, args.reverse_e_cutoff)?;
                let mmseqs_scores = decoy_scores::<BlastTable>(
                    &results,
                    MMSEQS_DECOY,
                    family,
                    args.reverse_e_cutoff,
                )?;
                let hmmer_scores = decoy_scores::<HmmerTable>(
                    &results,
                    HMMER_DECOY,
                    family,
                    args.reverse_e_cutoff,
                )?;

                let (Some(n), Some(m)) = (nail_scores, mmseqs_scores) else {
                    // a family without both tables tells us nothing comparative
                    skipped.fetch_add(1, Ordering::Relaxed);
                    return Ok(());
                };

                let mut line = format!("{family},{},{}", group("nail", &n), group("mmseqs", &m));
                if let Some(h) = hmmer_scores {
                    line.push(',');
                    line.push_str(&group("hmmer", &h));
                }

                writeln!(out.lock().expect("output mutex poisoned"), "{line}")?;
                Ok(())
            })
    })?;

    out.into_inner().expect("output mutex poisoned").flush()?;

    let skipped = skipped.load(Ordering::Relaxed);
    if skipped > 0 {
        eprintln!("skipped {skipped} families missing a nail or mmseqs table");
    }
    println!("wrote {}", out_path.display());
    Ok(())
}

/// The top decoy scores for one family and tool, plus how many decoys survived.
///
/// A reversed hit only counts as a decoy if the same (query, target) pair did
/// not also hit in the forward direction: reversal preserves composition, so a
/// genuine family member's reversal can score for reasons that are not chance.
fn decoy_scores<T>(
    results: &Path,
    run: &str,
    family: &str,
    e_cutoff: f64,
) -> anyhow::Result<Option<(Vec<f32>, usize)>>
where
    T: bioio::tbl::HitColumns,
{
    let fwd_path = results.join(format!("{run}.{family}.tbl"));
    let rev_path = results.join(format!("{run}.{family}.rev.tbl"));

    if !fwd_path.exists() || !rev_path.exists() {
        return Ok(None);
    }

    let mut fwd = HitTable::from_path::<_, T>(&fwd_path)
        .with_context(|| format!("failed to read {}", fwd_path.display()))?;
    let rev = HitTable::from_path::<_, T>(&rev_path)
        .with_context(|| format!("failed to read {}", rev_path.display()))?;

    fwd.hits.retain(|h| h.e_value <= e_cutoff);
    let real = fwd.to_map();

    let mut decoys: Vec<_> = rev
        .to_map()
        .into_iter()
        .filter(|(k, _)| !real.contains_key(k))
        .map(|(_, v)| v)
        .collect();

    decoys.sort_by(|a, b| {
        a.e_value
            .partial_cmp(&b.e_value)
            .expect("NaN in decoy e-values")
    });

    let scores: Vec<f32> = decoys
        .iter()
        .map(|h| h.score)
        // a family with fewer than N_SCORES decoys pads with zero, which
        // parse_cutoffs reads as "no usable cutoff"
        .chain(std::iter::repeat(0.0))
        .take(N_SCORES)
        .collect();

    Ok(Some((scores, decoys.len())))
}

/// One tool's group in a cutoffs line: `(tool,s1,...,sN,count)`.
fn group(tool: &str, (scores, count): &(Vec<f32>, usize)) -> String {
    let scores = scores
        .iter()
        .map(|s| format!("{s:.1}"))
        .collect::<Vec<_>>()
        .join(",");
    format!("({tool},{scores},{count})")
}

// --------------------------------------------------------------------- all

fn all(args: AllArgs) -> anyhow::Result<()> {
    reverse(ReverseArgs {
        place: args.place.clone(),
        shards: args.shards,
        threads: args.threads,
    })?;

    recruit(RecruitArgs {
        place: args.place.clone(),
        threads: args.threads,
        dry_run: false,
    })?;

    decoys(DecoysArgs {
        place: args.place.clone(),
        threads: args.threads,
    })?;

    search(SearchArgs {
        place: args.place.clone(),
        jobs: args.jobs,
    })?;

    learn(LearnArgs {
        place: args.place.clone(),
        reverse_e_cutoff: 1e-3,
        threads: args.threads,
    })
}

// ------------------------------------------------------------------- utils

/// A thread pool scoped to one phase, so it does not fight with any global one.
fn pool(threads: usize) -> anyhow::Result<rayon::ThreadPool> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .context("failed to build a thread pool")
}
