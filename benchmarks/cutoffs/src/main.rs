//! Learning per-family false-positive score cutoffs from reversed decoys.
//!
//! Reversing a protein keeps its composition but destroys its homology, so a
//! reversed sequence that still scores against a family is measuring noise.
//! Collecting those scores per family gives a threshold above which a hit is
//! unlikely to be chance.
//!
//! The set it reads is already reversed -- a `fixed` recipe under a
//! `reversed` tag, drawn the same way and written backwards -- so there is no
//! stage here that reverses anything.
//!
//! Doing that exhaustively would mean searching every family against every
//! reversed sequence. Instead it runs in two stages:
//!
//!   1. `recruit` — a cheap sweep of every family against the initial
//!      reversals, which finds the small subset that scores at all. What it
//!      finds are *recruits*.
//!   2. `reject` — an exhaustive pass, against a family's own recruits in both
//!      forms. A recruit whose *original* does not match the family is a
//!      *decoy*; one whose original does match is a *reject*, a reversed
//!      homolog rather than a piece of noise.
//!
//! `reject` is the stage the cutoffs are read from. It runs at high sensitivity
//! with the prefilter effectively disabled, so a decoy's score is its real
//! score rather than one truncated by `recruit`'s parameters. Searching the
//! originals is what separates the decoys from the rejects, and it matters
//! more than it looks: a reversal keeps a surprisingly high score against its
//! own original, so the recruits that look like the best decoys are exactly
//! the ones that may not be decoys at all.
//!
//! `gather` sits between them: it reads `recruit`'s hit tables and pulls the
//! sequences that hit out of the shard they came from. A record read there is
//! already the reversal, and reversing it again gives the original, so both
//! forms come out of one pass over the recruits.
//!
//! Both searching stages are `michi` pipelines, and they get there differently.
//! `recruit` is one big search per shard, so a shard's short chain unrolls
//! straight into steps. `search` is the opposite shape: many single-query
//! searches, each a chain of its own, and every tool here parallelises over
//! queries alone -- so a thread count above one buys nothing and the
//! parallelism has to be many families at once.
//!
//! A `Step` holds `Cmd`s rather than `Step`s, so a batch of ordered chains is
//! not something michi can be asked for. `search` transposes it: one step per
//! link of the chain, each batched across every family. Every family's link
//! still runs in order, since a step finishes before the next begins, and the
//! cost is a barrier per link rather than one straggler overall.
//!
//! `decoys` and `learn` run no tools, so they stay a plain rayon pool.

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
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
use libsail::tbl::{Hit, HitColumns, HitParser, Table};
use michi::{Cmd as PCmd, PipelineBuilder, Progress, Step, Table as PTable};
use util::tools::{hmmsearch, mmseqs, nail};
use util::{ledger, manifest, tbl};

use util::set::Set;

use util::cut;

/// Name of the calibration when none is given.
pub const DEFAULT_NAME: &str = "default";

/// What `gather` lays out, relative to its own root. These spell both the paths
/// the stages write and the cells the manifest carries.
//
// recruits, not decoys: a recruit only becomes a decoy once `reject` has shown
// that its original does not match the family that recruited it
const REVERSALS: &str = "reversals";
const ORIGINALS: &str = "originals";
const QUERIES: &str = "queries";

// recruitment only has to nominate candidates, so it runs a cheap sweep
// rather than the decoy stage's wide-open one
const RECRUIT_S: &str = "11.0";
const RECRUIT_MAX_SEQS: usize = 5000;

// wide open, so a decoy's score is its real score rather than one truncated
// by a prefilter. one thread each: the parallelism is in running many
// families at once, not in any one search
const DECOY_S: &str = "12.0";
const DECOY_MAX_SEQS: usize = 1_000_000_000;

// the run names both stages file their tables under. Each stage has its own
// results directory, so the stage does not need naming again in the file
const NAIL: &str = "nail";
const MMSEQS: &str = "mmseqs";
const HMMER: &str = "hmmer";

/// Which form of a family's recruits a search covered. The originals keep the
/// bare tool name; only the reversals need saying.
/// What the manifest calls the commands around a search, so a database build
/// is never charged to the tool that reads it.
//
// createdb and convert are not here: those commands come from
// `search`, which names its own stages
const DIRS: &str = "dirs";
const PROFILE: &str = "profile";
const CLEAN: &str = "clean";

/// The column that tells a search of the reversals from one of the originals.
//
// `form` rather than `direction`: nail's Forward algorithm owns that word here,
// and a reversal and its original are two forms of one sequence rather than
// two directions of anything
const FORM: &str = "form";

const ORIGINAL: &str = "orig";
const REVERSAL: &str = "rev";

/// The run name for one tool searching one form.
fn run_name(tool: &str, form: &str) -> String {
    match form {
        REVERSAL => format!("{tool}-{REVERSAL}"),
        _ => tool.to_string(),
    }
}

// ------------------------------------------------------------------ layout

/// Where every artifact of a calibration lives.
///
/// The calibration is three things wearing one name. The first two stages
/// build a set of decoys, which is a set like any other and goes where sets
/// go. The next searches it, which is a run. The last reads that run, which is
/// an analysis. Nothing here is a special kind of artifact; what makes it a
/// calibration rather than a benchmark is that its product is a data file to
/// be promoted by hand.
///
/// ```text
/// <set>/outputs/recruit/            recruit's tables
/// <set>/outputs/gather/
///   reversals/<family>.fa           the recruits, as they were hit
///   originals/<family>.fa           ... re-reversed, which is the original
///   queries/<family>/               the query set, per family
/// <set>/outputs/reject/             reject's tables
/// <set>/analysis/cutoffs/           what was learned
/// ```
struct Layout {
    /// What `gather` makes for the later stages to read.
    ///
    /// Outputs of this pipeline rather than a set of their own: a set is what
    /// `build-set` produces from sources under a recipe, and these come out of
    /// whatever `recruit` happened to score.
    gather: PathBuf,
    /// Where each of the two searching stages writes, and the scratch they
    /// share.
    recruit: PathBuf,
    reject: PathBuf,
    analysis: PathBuf,
    tmp: PathBuf,
    // resolved out of the source set's manifest once, here, so that the rest
    // of the calibration works in paths rather than in lookups -- and so a set
    // that cannot answer for one of them fails before any stage starts
    query_hmm: PathBuf,
    query_sto: PathBuf,
    query_db: PathBuf,
    targets: Vec<(String, PathBuf)>,
}

impl Layout {
    fn new(paths: &Paths) -> anyhow::Result<Self> {
        let set = Set::load_as(&paths.set, SHAPE)?;

        let first = set.units().next().context("the set names no units")?;
        let (query_hmm, query_sto, query_db) =
            (first.query_hmm()?, first.query_sto()?, first.query_db()?);

        let targets = set
            .units()
            .map(|u| Ok((u.name().to_string(), u.target()?)))
            .collect::<anyhow::Result<Vec<_>>>()?;

        Ok(Layout {
            gather: paths.gather.clone(),
            recruit: paths.recruit.clone(),
            reject: paths.reject.clone(),
            analysis: paths.analysis.clone(),
            tmp: paths.tmp.clone(),
            query_hmm,
            query_sto,
            query_db,
            targets,
        })
    }

    /// The shards the decoys are drawn from, as the source set named them.
    fn targets(&self) -> &[(String, PathBuf)] {
        &self.targets
    }

    fn query_hmm(&self) -> PathBuf {
        self.query_hmm.clone()
    }

    fn query_sto(&self) -> PathBuf {
        self.query_sto.clone()
    }

    /// The mmseqs profile db the build made from the query set. Reused here
    /// rather than rebuilt, since recruit searches the whole query set.
    fn query_db(&self) -> PathBuf {
        self.query_db.clone()
    }

    fn originals(&self) -> PathBuf {
        self.gather.join(ORIGINALS)
    }

    fn reversals(&self) -> PathBuf {
        self.gather.join(REVERSALS)
    }

    /// One directory per family, each holding a `query.hmm` and a `query.sto`
    /// -- the same shape as a ladder rung's query directory.
    fn queries(&self) -> PathBuf {
        self.gather.join(QUERIES)
    }

    fn recruit(&self) -> anyhow::Result<Stage> {
        self.stage("recruit")
    }

    fn reject(&self) -> anyhow::Result<Stage> {
        self.stage("search")
    }

    fn stage(&self, stage: &str) -> anyhow::Result<Stage> {
        let root = match stage {
            "recruit" => self.recruit.clone(),
            _ => self.reject.clone(),
        };

        Ok(Stage {
            root,
            tmp: self.tmp.join(stage),
        })
    }

    fn cutoffs_tbl(&self) -> anyhow::Result<PathBuf> {
        Ok(self.analysis.join("cutoffs.tbl"))
    }
}

/// One stage's output directory, in the shape every pipeline in this crate
/// writes: what ran, what it produced, and the scratch it wanted.
struct Stage {
    root: PathBuf,
    tmp: PathBuf,
}

impl Stage {
    fn results(&self) -> PathBuf {
        self.root.join("results")
    }

    fn tmp(&self) -> PathBuf {
        self.tmp.clone()
    }

    fn manifest(&self) -> PathBuf {
        self.root.join("manifest.tbl")
    }
}

// ---------------------------------------------------------------- commands

#[derive(Subcommand)]
enum Cmd {
    /// Search every family against the initial reversals to find recruits.
    Recruit(RecruitArgs),
    /// Lay out the recruits and their originals, and split the query set.
    Gather(GatherArgs),
    /// Search each family against its recruits, both forms, which splits them
    /// into decoys and rejects.
    Reject(RejectArgs),
    /// Turn the decoy scores into per-family cutoffs.
    Learn(LearnArgs),
    /// Run every stage in order.
    All(AllArgs),
}

/// Which label of paths.toml to work under.
#[derive(Parser, Debug, Clone)]
pub struct Where {
    /// Omit to list the labels
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,
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
pub struct GatherArgs {
    #[command(flatten)]
    pub place: Where,

    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Parser, Debug)]
pub struct RejectArgs {
    #[command(flatten)]
    pub place: Where,

    #[arg(long)]
    pub dry_run: bool,

    /// How many families to search at once. Each search is single-threaded,
    /// so this is the whole of the parallelism
    #[arg(short = 'j', long)]
    pub jobs: Option<usize>,
}

#[derive(Parser, Debug)]
pub struct LearnArgs {
    #[command(flatten)]
    pub place: Where,

    /// A recruit whose original hits at or below this E-value is a reject, and
    /// its reversal is left out of the decoy scores
    #[arg(short = 'e', default_value_t = 1e-3, value_name = "F")]
    pub reverse_e_cutoff: f64,

    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Parser, Debug)]
pub struct AllArgs {
    #[command(flatten)]
    pub place: Where,


    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,

    #[arg(short = 'j', long)]
    pub jobs: Option<usize>,
}

/// The set this calibrates against: the same deal a benchmark searches, drawn
/// backwards. Everything it makes out of that is its own output rather than a
/// second set.
pub const SHAPE: &util::set::Shape = &util::set::shape::REVERSED;

#[derive(Parser)]
#[command(name = "cutoffs", about = "per-family score cutoffs from reversed decoys")]
struct Cli {
    #[command(subcommand)]
    command: Cmd,
}

/// Where this calibration reads and writes, as one label of `paths.toml` names
/// them.
///
/// It is the one tool here that both reads a set and builds one, so it names
/// both: `set` is what it calibrates, `decoys` is what it makes out of that.
#[derive(serde::Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct Paths {
    pub set: PathBuf,
    pub gather: PathBuf,
    pub recruit: PathBuf,
    pub reject: PathBuf,
    pub analysis: PathBuf,
    pub tmp: PathBuf,
}

impl Paths {
    pub fn open(label: &str) -> anyhow::Result<Paths> {
        let file = util::paths::File::open(env!("CARGO_MANIFEST_DIR"))?;
        let p: Paths = file.get(label)?;

        Ok(Paths {
            set: file.at(p.set),
            gather: file.at(p.gather),
            recruit: file.at(p.recruit),
            reject: file.at(p.reject),
            analysis: file.at(p.analysis),
            tmp: file.at(p.tmp),
        })
    }

    pub fn listing(usage: &str) -> anyhow::Result<String> {
        Ok(util::paths::File::open(env!("CARGO_MANIFEST_DIR"))?.listing(usage))
    }
}

const USAGE: &str = "cutoffs <stage> --in <label>";

fn main() -> anyhow::Result<()> {
    let cmd = Cli::parse().command;

    let Some(label) = cmd.label() else {
        println!("{}", Paths::listing(USAGE)?);
        return Ok(());
    };
    let paths = Paths::open(label)?;

    run_cmd(cmd, &paths)
}

impl Cmd {
    fn label(&self) -> Option<&str> {
        let place = match self {
            Cmd::Recruit(a) => &a.place,
            Cmd::Gather(a) => &a.place,
            Cmd::Reject(a) => &a.place,
            Cmd::Learn(a) => &a.place,
            Cmd::All(a) => &a.place,
        };
        place.label.as_deref()
    }
}

fn run_cmd(cmd: Cmd, paths: &Paths) -> anyhow::Result<()> {
    match cmd {
        Cmd::Recruit(args) => recruit(args, paths),
        Cmd::Gather(args) => gather(args, paths),
        Cmd::Reject(args) => reject(args, paths),
        Cmd::Learn(args) => learn(args, paths),
        Cmd::All(args) => all(args, paths),
    }
}

// ----------------------------------------------------------------- recruit

fn recruit(args: RecruitArgs, paths: &Paths) -> anyhow::Result<()> {
    let layout = Layout::new(paths)?;

    let nail_bin = nail()?;
    let mmseqs_bin = mmseqs()?;

    let query_hmm = layout.query_hmm();
    let query_db = layout.query_db();

    let stage = layout.recruit()?;
    let (results, tmp) = (stage.results(), stage.tmp());

    let mut pl = PipelineBuilder::new().step(PCmd::new("mkdir").flag("-p").path(&results));

    for (idx, shard) in layout.targets() {
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
                    search::createdb(
                        &mmseqs_bin,
                        shard,
                        &target_db,
                        &idx.to_string(),
                        args.threads,
                    ),
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
                    .arg("--seed-mode", search::sweeps::SEED_MODE)
                    .arg("--mmseqs-max-seqs", RECRUIT_MAX_SEQS)
                    .arg("-E", search::EVALUE)
                    .arg(
                        "--tbl-out",
                        manifest::table_path(&results, NAIL, &idx.to_string()),
                    )
                    .flag("--allow-overwrite")
                    .path(&query_hmm)
                    .path(shard)])
                .name(format!("nail.{idx}")),
            )
            .step(
                Step::serial({
                    let cmds = search::Mmseqs {
                        bin: &mmseqs_bin,
                        query_db: &query_db,
                        target_db: &target_db,
                        aln_db: aln_db.clone(),
                        work: scratch.join("work"),
                        out: manifest::table_path(&results, MMSEQS, &idx.to_string()),
                        threads: args.threads,
                        s: Some(RECRUIT_S.to_string()),
                        max_seqs: Some(RECRUIT_MAX_SEQS),
                    }
                    .cmds();

                    [
                        cmds.search,
                        cmds.convert
                            .field(manifest::NAME, MMSEQS)
                            .field(manifest::TOOL, MMSEQS)
                            .field(manifest::SHARD, idx.to_string()),
                    ]
                })
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

fn gather(args: GatherArgs, paths: &Paths) -> anyhow::Result<()> {
    let layout = Layout::new(paths)?;
    let recruit_results = layout.recruit()?.results();

    let shard_list: Vec<String> = layout
        .targets()
        .iter()
        .map(|(name, _)| name.clone())
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
                    let tbl = Table::<HitParser<NailTable>>::open(&path)
                        .with_context(|| format!("failed to read {}", path.display()))?;
                    collect(&tbl, &mut map);

                    let path = manifest::table_path(&recruit_results, MMSEQS, shard);
                    let tbl = Table::<HitParser<BlastTable>>::open(&path)
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

    // ---- pull the recruits out of the shards, in both forms

    let decoy_dir = layout.originals();
    for dir in [&decoy_dir, &layout.reversals()] {
        if dir.exists() {
            std::fs::remove_dir_all(dir)?;
        }
    }
    std::fs::create_dir_all(&decoy_dir)?;

    let rev_dir = layout.reversals();
    std::fs::create_dir_all(&rev_dir)?;

    // one lock per family, over both directions: shards are read in parallel
    // and any of them may contribute to any family
    let handles: HashMap<&str, Mutex<(PathBuf, PathBuf)>> = families
        .iter()
        .map(|f| {
            let paths = (
                decoy_dir.join(format!("{f}.fa")),
                rev_dir.join(format!("{f}.fa")),
            );
            (f.as_str(), Mutex::new(paths))
        })
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

                let path = targets
                    .iter()
                    .find(|(name, _)| name == shard)
                    .map(|(_, path)| path.clone())
                    .with_context(|| format!("the set has no unit {shard:?}"))?;
                let shard_fa = IndexedFasta::open(&path)
                    .with_context(|| format!("failed to open {}", path.display()))?;

                // the shard is reversed, so a record read out of it is already
                // the reversal, and reversing it again is the original
                // one. both come out of the one pass, over the recruits alone
                // rather than over the shard
                let mut buffers: HashMap<&str, (Vec<u8>, Vec<u8>)> = HashMap::new();
                for mut rec in shard_fa.iter() {
                    let Some(fams) = by_target.get(rec.name_str()?).cloned() else {
                        continue;
                    };

                    let mut reversed = Vec::new();
                    rec.write_to(&mut reversed, DEFAULT_LINE_WIDTH)?;

                    rec.reverse();
                    let mut original = Vec::new();
                    rec.write_to(&mut original, DEFAULT_LINE_WIDTH)?;

                    for family in fams {
                        let (orig, rev) = buffers.entry(family).or_default();
                        orig.extend_from_slice(&original);
                        rev.extend_from_slice(&reversed);
                    }
                }

                for (family, (original, reversed)) in buffers {
                    let guard = handles
                        .get(family)
                        .with_context(|| format!("no handle for family {family}"))?
                        .lock()
                        .expect("family mutex poisoned");

                    let (orig_path, rev_path) = &*guard;
                    for (path, text) in [(orig_path, &original), (rev_path, &reversed)] {
                        let mut file = std::fs::OpenOptions::new()
                            .create(true)
                            .append(true)
                            .open(path)?;
                        file.write_all(text)?;
                    }
                }

                Ok(())
            })
    })?;

    // ---- split the query set, so each family can be searched on its own

    println!("splitting queries for {} families...", families.len());

    let queries = layout.queries();
    let hmm = cut::scatter_hmm(layout.query_hmm(), &families, &queries)?;
    let sto = cut::scatter_sto(layout.query_sto(), &families, &queries)?;

    if hmm != families.len() || sto != families.len() {
        bail!(
            "recruited {} families but found {hmm} in query.hmm and {sto} in query.sto",
            families.len()
        );
    }

    println!("wrote {}", layout.gather.display());
    Ok(())
}

/// Fold a hit table into a family to target-name map.
fn collect<C: HitColumns>(tbl: &Table<HitParser<C>>, map: &mut HashMap<String, Vec<String>>) {
    for hit in tbl.iter() {
        map.entry(hit.query.clone())
            .or_default()
            .push(hit.target.clone());
    }
}

// ------------------------------------------------------------------ search

fn reject(args: RejectArgs, paths: &Paths) -> anyhow::Result<()> {
    let layout = Layout::new(paths)?;
    let decoy_dir = layout.originals();

    if !decoy_dir.is_dir() {
        bail!(
            "no decoys in {}; run `cutoffs decoys` first",
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

    let rev_dir = layout.reversals();
    let queries = layout.queries();

    let stage = layout.reject()?;
    let results = stage.results();
    if results.exists() {
        std::fs::remove_dir_all(&results)?;
    }

    let tmp = stage.tmp();

    let nail_bin = nail()?;
    let mmseqs_bin = mmseqs()?;
    let hmmsearch_bin = hmmsearch()?;

    let jobs = args
        .jobs
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(4, |n| n.get()));

    println!("searching {} families, {jobs} at once...", families.len());

    // one command per family per link of the chain, and one Step per link.
    //
    // the searches here are single-query -- one family's profile against its
    // own decoys -- and every tool parallelises over queries, so a thread
    // count above 1 buys nothing and the parallelism has to come from running
    // many families at once. `batched` is that: `jobs` commands in flight,
    // each of them single-threaded.
    //
    // which is why the steps are links rather than families. a family's
    // commands have to run in order, and a Step holds Cmds rather than Steps,
    // so a batch of ordered chains is not a thing michi can be asked for. the
    // transpose is: every family's link k, then every family's link k+1. the
    // ordering each family needs still holds, since a step finishes before the
    // next one starts.
    let per_family = |f: &dyn Fn(&str) -> PCmd| -> Vec<PCmd> {
        families.iter().map(|family| f(family)).collect()
    };

    let scratch = |family: &str| tmp.join(family);
    let dir_scratch = |family: &str, direction: &str| scratch(family).join(direction);
    let query_db = |family: &str| scratch(family).join("queryDB");

    let mut pl = PipelineBuilder::new().step(
        Step::batched(
            jobs,
            per_family(&|family| {
                let mut cmd = PCmd::new("mkdir").name("dirs").flag("-p").path(&results);
                for direction in [ORIGINAL, REVERSAL] {
                    let d = dir_scratch(family, direction);
                    cmd = cmd.path(d.join("targetDB")).path(d.join("alnDB"));
                }
                cmd.field(manifest::STAGE, DIRS)
            }),
        )
        .name("dirs"),
    );

    // the query profile is the same for both directions, so it is built once
    // per family rather than once per search
    pl = pl
        .step(
            Step::batched(
                jobs,
                per_family(&|family| {
                    PCmd::new(&mmseqs_bin)
                        .name("convertmsa")
                        .sub("convertmsa")
                        .path(queries.join(family).join("query.sto"))
                        .path(scratch(family).join("msaDB"))
                        .arg("--identifier-field", 0)
                        .field(manifest::STAGE, PROFILE)
                }),
            )
            .name("convertmsa"),
        )
        .step(
            Step::batched(
                jobs,
                per_family(&|family| {
                    PCmd::new(&mmseqs_bin)
                        .name("msa2profile")
                        .sub("msa2profile")
                        .path(scratch(family).join("msaDB"))
                        .path(query_db(family))
                        .arg("--match-mode", 1)
                        .field(manifest::STAGE, PROFILE)
                }),
            )
            .name("msa2profile"),
        );

    for (direction, decoys) in [(ORIGINAL, &decoy_dir), (REVERSAL, &rev_dir)] {
        let target = |family: &str| decoys.join(format!("{family}.fa"));
        let hmm = |family: &str| queries.join(family).join("query.hmm");
        let table = |tool: &str, family: &str| {
            manifest::table_path(&results, &run_name(tool, direction), family)
        };

        // the family is the shard and the direction is the run, so a table is
        // named the way every other pipeline names one
        let run_of = |tool: &'static str, family: &str| {
            (
                run_name(tool, direction),
                tool,
                family.to_string(),
            )
        };

        pl = pl
            .step(
                Step::batched(
                    jobs,
                    per_family(&|family| {
                        let (name, tool, shard) = run_of(NAIL, family);
                        PCmd::new(&nail_bin)
                            .sub("search")
                            .arg("--mmseqs-path", &mmseqs_bin)
                            .arg("-t", 1)
                            .arg("--tmp-dir", dir_scratch(family, direction).join("nail"))
                            .arg("--mmseqs-s", DECOY_S)
                            .arg("--seed-mode", search::sweeps::SEED_MODE)
                            .arg("--mmseqs-max-seqs", DECOY_MAX_SEQS)
                            .arg("-E", search::EVALUE)
                            .flag("--allow-overwrite")
                            .arg("--tbl-out", table(NAIL, family))
                            .path(hmm(family))
                            .path(target(family))
                            .field(manifest::NAME, name)
                            .field(manifest::TOOL, tool)
                            .field(manifest::SHARD, shard)
                            .field(FORM, direction)
                    }),
                )
                .name(format!("nail.{direction}")),
            )
            .step(
                Step::batched(
                    jobs,
                    per_family(&|family| {
                        search::createdb(
                            &mmseqs_bin,
                            &target(family),
                            &dir_scratch(family, direction).join("targetDB/targetDB"),
                            family,
                            1,
                        )
                    }),
                )
                .name(format!("createdb.{direction}")),
            )
            .step(
                Step::batched(
                    jobs,
                    per_family(&|family| {
                        let (name, tool, shard) = run_of(MMSEQS, family);
                        let d = dir_scratch(family, direction);
                        search::Mmseqs {
                            bin: &mmseqs_bin,
                            query_db: &query_db(family),
                            target_db: &d.join("targetDB/targetDB"),
                            aln_db: d.join("alnDB/alnDB"),
                            work: d.join("work"),
                            out: table(MMSEQS, family),
                            threads: 1,
                            s: Some(DECOY_S.to_string()),
                            max_seqs: Some(DECOY_MAX_SEQS),
                        }
                        .cmds()
                        .search
                        .field(manifest::NAME, name)
                        .field(manifest::TOOL, tool)
                        .field(manifest::SHARD, shard)
                        .field(FORM, direction)
                    }),
                )
                .name(format!("mmseqs.{direction}")),
            )
            .step(
                Step::batched(
                    jobs,
                    per_family(&|family| {
                        let d = dir_scratch(family, direction);
                        search::Mmseqs {
                            bin: &mmseqs_bin,
                            query_db: &query_db(family),
                            target_db: &d.join("targetDB/targetDB"),
                            aln_db: d.join("alnDB/alnDB"),
                            work: d.join("work"),
                            out: table(MMSEQS, family),
                            threads: 1,
                            s: Some(DECOY_S.to_string()),
                            max_seqs: Some(DECOY_MAX_SEQS),
                        }
                        .cmds()
                        .convert
                        .field(manifest::SHARD, family)
                    }),
                )
                .name(format!("convert.{direction}")),
            )
            .step(
                Step::batched(
                    jobs,
                    // hmmsearch is the one command here not built through
                    // `search`. search::hmmer returns a whole Step, batched
                    // over the parts of one split query and followed by a cat
                    // that joins their tables. this searches one family per
                    // command and batches over families instead, and each
                    // family's table is its own, so there is no split to cut
                    // and nothing to concatenate
                    per_family(&|family| {
                        let (name, tool, shard) = run_of(HMMER, family);
                        PCmd::new(&hmmsearch_bin)
                            .arg("--cpu", 1)
                            .arg("-E", search::EVALUE)
                            .arg("-o", "/dev/null")
                            .arg("--tblout", table(HMMER, family))
                            .arg(
                                "--domtblout",
                                manifest::dom_path(&results, &run_name(HMMER, direction), family),
                            )
                            .path(hmm(family))
                            .path(target(family))
                            .field(manifest::NAME, name)
                            .field(manifest::TOOL, tool)
                            .field(manifest::SHARD, shard)
                            .field(FORM, direction)
                    }),
                )
                .name(format!("hmmer.{direction}")),
            );
    }

    pl = pl.step(
        Step::batched(
            jobs,
            per_family(&|family| {
                PCmd::new("rm")
                    .name("clean")
                    .flag("-rf")
                    .path(scratch(family))
                    .field(manifest::STAGE, CLEAN)
            }),
        )
        .name("clean"),
    );

    let pipeline = pl
        .stderr_dir(tmp.join("stderr"))
        .sink(Progress::new())
        .sink(PTable::new(stage.manifest()))
        .build()
        .context("failed to build the search")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    ledger::clear(&stage.root);
    pipeline.run()?;
    ledger::record(&stage.root)
}

// ------------------------------------------------------------------- learn

/// How many decoy scores are recorded per family per tool.
const N_SCORES: usize = 5;

/// The tools a calibration scores, in the order `cutoffs.tbl` writes them.
const TOOLS: [&str; 3] = [NAIL, MMSEQS, HMMER];

fn learn(args: LearnArgs, paths: &Paths) -> anyhow::Result<()> {
    let layout = Layout::new(paths)?;
    let stage = layout.reject()?;
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
                let nail =
                    decoy_scores::<NailTable>(&results, NAIL, family, args.reverse_e_cutoff)?;
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

    let out_path = layout.cutoffs_tbl()?;
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
/// whose original does not also hit: reversal preserves composition, so a
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
    let table = |direction| manifest::table_path(results, &run_name(tool, direction), family);

    let (orig_path, rev_path) = (table(ORIGINAL), table(REVERSAL));

    if !orig_path.exists() || !rev_path.exists() {
        return Ok(None);
    }

    let orig = Table::<HitParser<T>>::open(&orig_path)
        .with_context(|| format!("failed to read {}", orig_path.display()))?;
    let rev = Table::<HitParser<T>>::open(&rev_path)
        .with_context(|| format!("failed to read {}", rev_path.display()))?;

    let real: HashSet<(&str, &str)> = orig
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

fn all(args: AllArgs, paths: &Paths) -> anyhow::Result<()> {
    recruit(RecruitArgs {
        place: args.place.clone(),
        threads: args.threads,
        dry_run: false,
    }, paths)?;

    gather(GatherArgs {
        place: args.place.clone(),
        threads: args.threads,
    }, paths)?;

    reject(
        RejectArgs {
            place: args.place.clone(),
            dry_run: false,
            jobs: args.jobs,
        },
        paths,
    )?;

    learn(
        LearnArgs {
            place: args.place.clone(),
            reverse_e_cutoff: 1e-3,
            threads: args.threads,
        },
        paths,
    )
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
