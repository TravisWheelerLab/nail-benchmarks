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
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{Context, bail};
use clap::{Parser, Subcommand};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use libsail::collection::Iterable;
use libsail::seq::fasta::{DEFAULT_LINE_WIDTH, IndexedFasta};
use libsail::tbl::blast::BlastTable;
use libsail::tbl::hmmer::HmmerTable;
use libsail::tbl::nail::NailTable;
use libsail::tbl::{Hit, HitColumns, Table};
use michi::{Cmd as PCmd, PipelineBuilder, Progress, Step, Table as PTable};
use util::tools::{hmmsearch, mmseqs, nail};
use util::{manifest, tbl};

use crate::cut;
use crate::inputs::{self, shards};

/// Name of the calibration when none is given.
pub const DEFAULT_NAME: &str = "default";

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

// the run names both stages file their tables under. Each stage has its own
// results directory, so the stage does not need naming again in the file
const NAIL: &str = "nail";
const MMSEQS: &str = "mmseqs";
const HMMER: &str = "hmmer";

/// Which direction of a family's decoys a search covered. The forward run keeps
/// the bare tool name; only the reversed one needs saying.
const FORWARD: &str = "fwd";
const REVERSE: &str = "rev";

/// The run name for one tool searching one direction.
fn run_name(tool: &str, direction: &str) -> String {
    match direction {
        REVERSE => format!("{tool}-{REVERSE}"),
        _ => tool.to_string(),
    }
}

// ------------------------------------------------------------------ layout

/// Where every artifact of a calibration lives.
///
/// The whole thing hangs off one directory at the crate root rather than under
/// `outputs/`, so several calibrations — different shard counts, different
/// parameters — can sit side by side and be deleted as a unit. What it produces
/// is a data file to be promoted by hand, not a benchmark's results.
///
/// Inside, it has the shape the rest of the crate has: what a stage generates
/// for a later stage to read is under `inputs/`, and what a stage produces for
/// its own sake is under `outputs/`.
///
/// ```text
/// cutoffs/<name>/inputs/targets-rev/<i>.fa        the reversed shards
/// cutoffs/<name>/inputs/decoys/<family>.fa        what recruited, un-reversed
/// cutoffs/<name>/inputs/decoys-rev/<family>.fa    ... and reversed again
/// cutoffs/<name>/inputs/queries/<family>/         the query set, per family
/// cutoffs/<name>/outputs/recruit/                 stage 1's tables
/// cutoffs/<name>/outputs/search/                  stage 2's tables
/// cutoffs/<name>/outputs/cutoffs.tbl              what was learned
/// ```
///
/// The inputs it starts from are the fixed set every recall pipeline reads.
struct Layout {
    root: PathBuf,
}

impl Layout {
    fn new(name: &str) -> anyhow::Result<Self> {
        let set = inputs::fixed::dir();
        if !set.is_dir() {
            bail!(
                "{} does not exist; run `mgy build fixed` first",
                set.display()
            );
        }

        // a name, not a path: joining an unchecked one would let a calibration
        // be written anywhere on the machine
        if name.is_empty()
            || name == ".."
            || name.contains(std::path::MAIN_SEPARATOR)
            || name.contains('/')
        {
            bail!("--name {name:?} is a directory name, not a path");
        }

        Ok(Layout {
            root: crate::dir().join("cutoffs").join(name),
        })
    }

    fn inputs(&self) -> PathBuf {
        self.root.join("inputs")
    }

    fn outputs(&self) -> PathBuf {
        self.root.join("outputs")
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
        self.inputs().join("targets-rev")
    }

    fn decoys(&self) -> PathBuf {
        self.inputs().join("decoys")
    }

    fn decoys_rev(&self) -> PathBuf {
        self.inputs().join("decoys-rev")
    }

    /// One directory per family, each holding a `query.hmm` and a `query.sto`
    /// -- the same shape as a ladder rung's query directory.
    fn queries(&self) -> PathBuf {
        self.inputs().join("queries")
    }

    fn recruit(&self) -> Stage {
        Stage {
            root: self.outputs().join("recruit"),
        }
    }

    fn search(&self) -> Stage {
        Stage {
            root: self.outputs().join("search"),
        }
    }

    fn cutoffs_tbl(&self) -> PathBuf {
        self.outputs().join("cutoffs.tbl")
    }
}

/// One stage's output directory, in the shape every pipeline in this crate
/// writes: what ran, what it produced, and the scratch it wanted.
struct Stage {
    root: PathBuf,
}

impl Stage {
    fn results(&self) -> PathBuf {
        self.root.join("results")
    }

    fn tmp(&self) -> PathBuf {
        self.root.join("tmp")
    }

    fn manifest(&self) -> PathBuf {
        self.root.join("manifest.tbl")
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

/// Which calibration to work in.
#[derive(Parser, Debug, Clone)]
pub struct Where {
    /// Names the calibration, under benchmarks/mgy/cutoffs/
    #[arg(long, default_value = DEFAULT_NAME)]
    pub name: String,
}

#[derive(Parser, Debug)]
pub struct ReverseArgs {
    #[command(flatten)]
    pub place: Where,

    /// Reverse only the first N shards. This is the one place the size of the
    /// calibration set is decided; every later stage uses whatever is here
    #[arg(short = 'n', long)]
    pub shards: Option<usize>,

    /// Threads for the reversal itself
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Parser, Debug)]
pub struct RecruitArgs {
    #[command(flatten)]
    pub place: Where,

    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,

    /// List the commands and exit without executing anything
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
    /// so this is the whole of the parallelism
    #[arg(short = 'j', long)]
    pub jobs: Option<usize>,
}

#[derive(Parser, Debug)]
pub struct LearnArgs {
    #[command(flatten)]
    pub place: Where,

    /// Forward hits at or below this E-value are treated as real, and their
    /// reversed counterparts are excluded from the decoy scores
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

/// Write a copy of `src` with every sequence reversed. Reversed sequences keep
/// the composition of the original but destroy its homology, which makes them
/// usable as decoys when calibrating score cutoffs.
fn write_reversed(src: &Path, dst: &Path) -> anyhow::Result<()> {
    if let Some(dir) = dst.parent() {
        std::fs::create_dir_all(dir)?;
    }

    let fa = IndexedFasta::open(src)
        .with_context(|| format!("failed to open {}", src.display()))?;
    let mut out = BufWriter::new(
        std::fs::File::create(dst)
            .with_context(|| format!("failed to create {}", dst.display()))?,
    );

    for mut rec in fa.iter() {
        rec.reverse();
        rec.write_to(&mut out, DEFAULT_LINE_WIDTH)?;
    }

    out.flush()?;
    Ok(())
}

fn reverse(args: ReverseArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.name)?;
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
                write_reversed(path, &dst.join(format!("{i}.fa")))
                    .with_context(|| format!("failed to reverse shard {i}"))
            })
    })?;

    println!("reversed {} shards", found.len());
    Ok(())
}

// ----------------------------------------------------------------- recruit

fn recruit(args: RecruitArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.name)?;
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

    let stage = layout.recruit();
    let (results, tmp) = (stage.results(), stage.tmp());

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
                    .field(manifest::NAME, NAIL)
                    .field(manifest::TOOL, NAIL)
                    .field(manifest::SHARD, idx.to_string())
                    .arg("--mmseqs-path", &mmseqs_bin)
                    .arg("-t", args.threads)
                    .arg("--tmp-dir", scratch.join("nail"))
                    .arg("--mmseqs-s", RECRUIT_S)
                    .arg("--mmseqs-max-seqs", RECRUIT_MAX_SEQS)
                    .arg("-E", EVALUE)
                    .arg(
                        "--tbl-out",
                        manifest::table_path(&results, NAIL, &idx.to_string()),
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
                        .field(manifest::NAME, MMSEQS)
                        .field(manifest::TOOL, MMSEQS)
                        .field(manifest::SHARD, idx.to_string())
                        .arg("--format-mode", 0)
                        .path(&query_db)
                        .path(&target_db)
                        .path(&aln_db)
                        .path(manifest::table_path(&results, MMSEQS, &idx.to_string())),
                ])
                .name(format!("mmseqs.{idx}")),
            )
            .step(PCmd::new("rm").name("clean").flag("-rf").path(&scratch));
    }

    let pipeline = pl
        .stderr_dir(tmp.join("stderr"))
        .sink(Progress::new())
        .sink(PTable::new(stage.manifest()))
        .build()?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    pipeline.run()
}

// ------------------------------------------------------------------ decoys

fn decoys(args: DecoysArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.name)?;
    let recruit_results = layout.recruit().results();

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

                    let path = manifest::table_path(&recruit_results, NAIL, shard);
                    let tbl = Table::<NailTable>::open(&path)
                        .with_context(|| format!("failed to read {}", path.display()))?;
                    collect(&tbl, &mut map);

                    let path = manifest::table_path(&recruit_results, MMSEQS, shard);
                    let tbl = Table::<BlastTable>::open(&path)
                        .with_context(|| format!("failed to read {}", path.display()))?;
                    collect(&tbl, &mut map);

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
                let shard_fa = IndexedFasta::open(&path)
                    .with_context(|| format!("failed to open {}", path.display()))?;

                let mut buffers: HashMap<&str, Vec<u8>> = HashMap::new();
                for rec in shard_fa.iter() {
                    let Some(fams) = by_target.get(rec.name_str()?) else {
                        continue;
                    };
                    for family in fams {
                        let buf = buffers.entry(family).or_default();
                        rec.write_to(buf, DEFAULT_LINE_WIDTH)?;
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
                    file.write_all(&text)?;
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
            // plain `<family>.fa`, not `<family>.rev.fa`: the directory already
            // says these are reversed, the same call `reverse` makes for shards
            write_reversed(
                &decoy_dir.join(format!("{f}.fa")),
                &rev_dir.join(format!("{f}.fa")),
            )
            .with_context(|| format!("failed to reverse decoys for {f}"))
        })
    })?;

    // ---- split the query set, so each family can be searched on its own

    println!("splitting queries for {} families...", families.len());

    let queries = layout.queries();
    let hmm = cut::explode_hmm(layout.query_hmm(), &families, &queries)?;
    let sto = cut::explode_sto(layout.query_sto(), &families, &queries)?;

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
fn collect<C: HitColumns>(tbl: &Table<C>, map: &mut HashMap<String, Vec<String>>) {
    for hit in tbl.iter() {
        map.entry(hit.query.clone())
            .or_default()
            .push(hit.target.clone());
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
    let layout = Layout::new(&args.place.name)?;
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
    let queries = layout.queries();

    let stage = layout.search();
    let results = stage.results();
    if results.exists() {
        std::fs::remove_dir_all(&results)?;
    }
    std::fs::create_dir_all(&results)?;

    let tmp = stage.tmp();
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

                let hmm = queries.join(family).join("query.hmm");
                let sto = queries.join(family).join("query.sto");

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

                // the family is the shard and the direction is the run, so a
                // table is named the way every other pipeline names one
                let directions = [
                    (FORWARD, decoy_dir.join(format!("{family}.fa"))),
                    (REVERSE, rev_dir.join(format!("{family}.fa"))),
                ];

                for (direction, target) in directions {
                    let dir_scratch = scratch.join(direction);
                    std::fs::create_dir_all(&dir_scratch)?;

                    let table = |tool: &str| {
                        manifest::table_path(&results, &run_name(tool, direction), family)
                    };

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
                            .arg(table(NAIL))
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
                            .arg(table(MMSEQS))
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
                            .arg(table(HMMER))
                            .arg("--domtblout")
                            .arg(manifest::dom_path(
                                &results,
                                &run_name(HMMER, direction),
                                family,
                            ))
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

    write_manifest(&stage.manifest(), &families)?;

    println!("wrote {}", results.display());
    Ok(())
}

/// Record what this stage searched, in the shape `parse` reads elsewhere.
///
/// Assembled by hand rather than by a [`PTable`] sink: the searches run on a
/// rayon pool and not through a pipeline. Only the rows a reader needs are
/// here -- a name, what produced it, which family, which direction, and that it
/// finished. Every command is checked by `run`, so reaching this point means
/// they all did.
fn write_manifest(path: &Path, families: &[String]) -> anyhow::Result<()> {
    let headers = ["name", "tool", "shard", "direction", "exit"]
        .map(str::to_string)
        .to_vec();

    let rows: Vec<Vec<String>> = families
        .iter()
        .flat_map(|family| {
            [FORWARD, REVERSE].into_iter().flat_map(move |direction| {
                [NAIL, MMSEQS, HMMER].into_iter().map(move |tool| {
                    vec![
                        run_name(tool, direction),
                        tool.to_string(),
                        family.clone(),
                        direction.to_string(),
                        "0".to_string(),
                    ]
                })
            })
        })
        .collect();

    tbl::write(
        path,
        tbl::Table {
            meta: "",
            headers: &headers,
            rows: &rows,
            ragged_last: false,
        },
    )
}

// ------------------------------------------------------------------- learn

/// How many decoy scores are recorded per family per tool.
const N_SCORES: usize = 5;

/// The tools a calibration scores, in the order `cutoffs.tbl` writes them.
const TOOLS: [&str; 3] = [NAIL, MMSEQS, HMMER];

fn learn(args: LearnArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.name)?;
    let stage = layout.search();
    let results = stage.results();

    // what the search stage actually did, rather than what is on disk: a run
    // that died leaves a table behind, and a half-written one reads as a family
    // with fewer decoys than it has
    let searched = manifest::Manifest::read(&stage.manifest())?;

    let mut families: Vec<String> = searched
        .runs()
        .filter_map(|row| row.get(manifest::SHARD).map(str::to_string))
        .collect();
    families.sort_unstable();
    families.dedup();

    let failed = searched.failed().count();
    if failed > 0 {
        bail!(
            "{failed} searches in {} did not finish; re-run `mgy cutoffs search`",
            stage.manifest().display()
        );
    }

    if families.is_empty() {
        bail!("no finished searches in {}", stage.manifest().display());
    }

    let skipped = AtomicUsize::new(0);

    // collected rather than written as they finish, so the file is in family
    // order however the pool interleaves
    let pool = pool(args.threads)?;
    let rows: Vec<Option<Vec<String>>> = pool.install(|| {
        families
            .par_iter()
            .map(|family| -> anyhow::Result<Option<Vec<String>>> {
                let nail = decoy_scores::<NailTable>(&results, NAIL, family, args.reverse_e_cutoff)?;
                let mmseqs =
                    decoy_scores::<BlastTable>(&results, MMSEQS, family, args.reverse_e_cutoff)?;
                let hmmer =
                    decoy_scores::<HmmerTable>(&results, HMMER, family, args.reverse_e_cutoff)?;

                if nail.is_none() || mmseqs.is_none() {
                    // a family without both tables tells us nothing comparative
                    skipped.fetch_add(1, Ordering::Relaxed);
                    return Ok(None);
                }

                let mut row = vec![family.clone()];
                for scores in [nail, mmseqs, hmmer] {
                    row.extend(cells(scores.as_ref()));
                }

                Ok(Some(row))
            })
            .collect::<anyhow::Result<Vec<_>>>()
    })?;

    let rows: Vec<Vec<String>> = rows.into_iter().flatten().collect();

    let mut headers = vec!["family".to_string()];
    for tool in TOOLS {
        headers.extend((1..=N_SCORES).map(|i| format!("{tool}_{i}")));
        headers.push(format!("{tool}_n"));
    }

    let out_path = layout.cutoffs_tbl();
    tbl::write(
        &out_path,
        tbl::Table {
            meta: "",
            headers: &headers,
            rows: &rows,
            ragged_last: false,
        },
    )?;

    let skipped = skipped.load(Ordering::Relaxed);
    if skipped > 0 {
        eprintln!("skipped {skipped} families missing a nail or mmseqs table");
    }
    println!("wrote {}", out_path.display());
    Ok(())
}

/// One tool's cells: its `N_SCORES` decoy scores and how many decoys there
/// were. A tool that scored none writes zeros, which is what a reader already
/// treats as "no usable cutoff".
fn cells(scores: Option<&(Vec<f32>, usize)>) -> Vec<String> {
    match scores {
        Some((scores, count)) => scores
            .iter()
            .map(|s| format!("{s:.1}"))
            .chain(std::iter::once(count.to_string()))
            .collect(),
        None => std::iter::repeat_n("0.0".to_string(), N_SCORES)
            .chain(std::iter::once("0".to_string()))
            .collect(),
    }
}

/// The top decoy scores for one family and tool, plus how many decoys survived.
///
/// A reversed hit only counts as a decoy if the same (query, target) pair did
/// not also hit in the forward direction: reversal preserves composition, so a
/// genuine family member's reversal can score for reasons that are not chance.
fn decoy_scores<T>(
    results: &Path,
    tool: &str,
    family: &str,
    e_cutoff: f64,
) -> anyhow::Result<Option<(Vec<f32>, usize)>>
where
    T: HitColumns,
{
    let table = |direction| {
        manifest::table_path(results, &run_name(tool, direction), family)
    };

    let (fwd_path, rev_path) = (table(FORWARD), table(REVERSE));

    if !fwd_path.exists() || !rev_path.exists() {
        return Ok(None);
    }

    let fwd = Table::<T>::open(&fwd_path)
        .with_context(|| format!("failed to read {}", fwd_path.display()))?;
    let rev = Table::<T>::open(&rev_path)
        .with_context(|| format!("failed to read {}", rev_path.display()))?;

    let real: HashSet<(&str, &str)> = fwd
        .iter()
        .filter(|h| h.e_value <= e_cutoff)
        .map(|h| (h.query.as_str(), h.target.as_str()))
        .collect();

    // one entry per pair, the last row for it winning, so a pair reported
    // twice counts as one decoy
    let mut by_pair: HashMap<(&str, &str), &Hit> = HashMap::new();
    for hit in rev.iter() {
        by_pair.insert((hit.query.as_str(), hit.target.as_str()), hit);
    }

    let mut decoys: Vec<&Hit> = by_pair
        .into_iter()
        .filter(|(k, _)| !real.contains(k))
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
        // `cutoffs` in scores.rs reads as "no usable cutoff"
        .chain(std::iter::repeat(0.0))
        .take(N_SCORES)
        .collect();

    Ok(Some((scores, decoys.len())))
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

/// A thread pool scoped to one phase, so it never shares threads with a global
/// one.
fn pool(threads: usize) -> anyhow::Result<rayon::ThreadPool> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .context("failed to build a thread pool")
}
