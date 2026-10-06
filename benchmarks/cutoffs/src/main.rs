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
//! reversed sequence. Instead it runs in stages:
//!
//!   1. `recruit` — a cheap sweep of nail and mmseqs over the initial
//!      reversals, which finds the small subset that scores at all. What it
//!      finds are *recruits*.
//!   2. `gather` — pulls each recruit out of the shard it came from. A record
//!      read there is already the reversal, and reversing it again gives the
//!      *original*, so both forms come out of one pass.
//!   3. `reject` — one hmmsearch of each family against its recruits'
//!      originals, at the `-Z` of one shard. A recruit whose original the
//!      family hits is a *reject*, a reversed homolog rather than a piece of
//!      noise; the rest are *decoys*.
//!   4. `fill` — nail, mmseqs and hmmsearch over the decoys, with the
//!      prefilter effectively off and `-E 1e9`, so that every tool has a
//!      score for every decoy and not only for the ones it recruited. It
//!      rescores nothing: a tool's score for a pair does not move with these
//!      settings. What the pass adds is the pairs a tool never reached in
//!      `recruit`.
//!   5. `learn` — the top scores per family per tool, into `cutoffs.tbl`.
//!
//! Rejection matters more than it looks: a reversal keeps a surprisingly high
//! score against its own original, so the recruits that look like the best
//! decoys are exactly the ones that may not be decoys at all. The judge
//! searches at one shard's `-Z` so that a reject means what it would in the
//! search a cutoff is applied to rather than in the small file of one
//! family's recruits.
//!
//! The searching stages are `michi` pipelines, and they get there differently.
//! `recruit` is one big search per shard, so a shard's short chain unrolls
//! straight into steps. `reject` and `fill` under `fanout` are the opposite
//! shape: many single-query searches, each a chain of its own, and every tool
//! here parallelises over queries alone -- so a thread count above one buys
//! nothing and the parallelism has to be many families at once.
//!
//! A `Step` holds `Cmd`s rather than `Step`s, so a batch of ordered chains is
//! not something michi can be asked for. `fill` transposes it: one step per
//! link of the chain, each batched across every family. Every family's link
//! still runs in order, since a step finishes before the next begins, and the
//! cost is a barrier per link rather than one straggler overall.
//!
//! `gather` and `learn` run no tools, so they stay a plain rayon pool.

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{Context, bail};
use clap::{Parser, Subcommand};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use libsail::collection::Iterable;
use libsail::format::Format;
use libsail::index::Reader;
use libsail::seq::fasta::{DEFAULT_LINE_WIDTH, IndexedFasta};
use libsail::tbl::blast::BlastTable;
use libsail::tbl::hmmer::HmmerTable;
use libsail::tbl::nail::NailTable;
use libsail::tbl::{Hit, HitColumns, HitParser, Table};
use michi::{Cmd as PCmd, PipelineBuilder, Progress, Step, Table as PTable};
use util::tools::{hmmsearch, mmseqs, nail};
use util::{ledger, manifest};

use util::set::Set;

use util::cut;

/// How `reject` and `fill` search: one process per family, or one per tool.
//
// the two are meant to agree exactly: same searches, same pairs, same cutoffs.
// they differ only in how many processes it takes to get there, and which is
// faster is an open question -- see the pinned note in CLAUDE.md
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Strategy {
    /// One search per family: its profile against its own recruits. Cheap per
    /// search, and there is one for every family in the set.
    Fanout,
    /// One search per tool: every family against every recruit at once, with
    /// each hit kept only for the family whose sequence it is.
    Union,
}

/// What `gather` lays out, relative to its own root. These spell both the paths
/// the stages write and the cells the manifest carries.
//
// recruits, not decoys: a recruit only becomes a decoy once `reject` has shown
// that its original does not match the family that recruited it
const REVERSALS: &str = "reversals";
const ORIGINALS: &str = "originals";
const QUERIES: &str = "queries";

/// The shard a union search files its one table under.
const ALL: &str = "all";

// recruitment only has to nominate candidates, so it runs a cheap sweep
// rather than the decoy stage's wide-open one
const RECRUIT_S: &str = "11.0";
const RECRUIT_MAX_SEQS: usize = 5000;

// wide open, so every decoy a tool can align gets a row: the point of the
// fill is a score from every tool for every decoy, including the ones that
// tool never recruited, and a reporting threshold would hide the weak ones.
// one E-value for the three tools rather than a -Z for the two that take one
const FILL_S: &str = "12.0";
const FILL_MAX_SEQS: usize = 1_000_000_000;
const FILL_E: &str = "1e9";

/// The E-value at which the judge's hit on a recruit's original makes it a
/// reject.
const REJECT_E: f64 = 1e-3;

// the run names both stages file their tables under. Each stage has its own
// results directory, so the stage does not need naming again in the file
const NAIL: &str = "nail";
const MMSEQS: &str = "mmseqs";
const HMMER: &str = "hmmer";

/// The run that searches the originals: hmmsearch, judging every tool's
/// recruits.
const JUDGE: &str = "judge";

/// Where `fill` lays out what it searched: the recruits the judge kept.
const DECOYS: &str = "decoys";

/// What the manifest calls the commands around a search, so a database build
/// is never charged to the tool that reads it.
//
// createdb and convert are not here: those commands come from
// `util::search`, which names its own stages
const DIRS: &str = "dirs";
const PROFILE: &str = "profile";
const CLEAN: &str = "clean";

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
    /// Where each of the three searching stages writes, and the scratch they
    /// share.
    recruit: PathBuf,
    reject: PathBuf,
    fill: PathBuf,
    analysis: PathBuf,
    tmp: PathBuf,
    // resolved out of the source set's manifest once, here, so that the rest
    // of the calibration works in paths rather than in lookups -- and so a set
    // that cannot answer for one of them fails before any stage starts
    query_hmm: PathBuf,
    query_sto: PathBuf,
    query_db: PathBuf,
    targets: Vec<(String, PathBuf)>,
    /// The sequences in one shard, which is the `-Z` the judge searches at.
    //
    // a cutoff learned here is applied to a search of one shard this size,
    // so a reject has to mean what it would mean there rather than in the
    // small file of one family's recruits. every shard of a fixed deal holds
    // the same count, and `new` refuses a set where they differ rather than
    // averaging them
    seqs: u64,
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

        let seqs = first.number("seqs")?;
        for unit in set.units() {
            let n = unit.number("seqs")?;
            if n != seqs {
                bail!(
                    "the set at {} has shards of different sizes ({seqs} and {n} sequences), so there is no one -Z for the judge",
                    paths.set.display()
                );
            }
        }

        Ok(Layout {
            gather: paths.gather.clone(),
            recruit: paths.recruit.clone(),
            reject: paths.reject.clone(),
            fill: paths.fill.clone(),
            analysis: paths.analysis.clone(),
            tmp: paths.tmp.clone(),
            query_hmm,
            query_sto,
            query_db,
            targets,
            seqs,
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

    /// One directory per family, each holding a `query.hmm` and a `query.sto`.
    fn queries(&self) -> PathBuf {
        self.gather.join(QUERIES)
    }

    fn recruit(&self) -> anyhow::Result<Stage> {
        self.stage("recruit")
    }

    fn reject(&self) -> anyhow::Result<Stage> {
        self.stage("reject")
    }

    fn fill(&self) -> anyhow::Result<Stage> {
        self.stage("fill")
    }

    /// The (family, sequence) pairs the judge rejected, written by `reject`.
    fn rejects_tbl(&self) -> PathBuf {
        self.reject.join("rejects.tbl")
    }

    /// One fasta of decoys per family, written by `fill`.
    fn decoys(&self) -> PathBuf {
        self.fill.join(DECOYS)
    }

    /// The (family, sequence) pairs `fill` searched, written by `fill`.
    fn decoys_tbl(&self) -> PathBuf {
        self.fill.join("decoys.tbl")
    }

    fn stage(&self, stage: &str) -> anyhow::Result<Stage> {
        let root = match stage {
            "recruit" => self.recruit.clone(),
            "fill" => self.fill.clone(),
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
    /// Search each family against its recruits' originals with hmmsearch,
    /// which splits them into decoys and rejects.
    Reject(RejectArgs),
    /// Search every tool against the decoys, so each has a score for every one.
    Fill(FillArgs),
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

    /// Threads per search
    #[arg(short, long)]
    pub threads: Option<usize>,

    /// List the commands and exit without executing anything
    #[arg(long)]
    pub dry_run: bool,
}

#[derive(Parser, Debug)]
pub struct GatherArgs {
    #[command(flatten)]
    pub place: Where,

    /// What `reject` and `fill` will run with. `union` searches the whole
    /// query set at once, so the per-family split is not laid out for it
    #[arg(long, value_enum, default_value_t = Strategy::Fanout)]
    pub strategy: Strategy,

    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Parser, Debug)]
pub struct RejectArgs {
    #[command(flatten)]
    pub place: Where,

    /// How to search. The two agree; they cost differently
    #[arg(long, value_enum, default_value_t = Strategy::Fanout)]
    pub strategy: Strategy,

    #[arg(long)]
    pub dry_run: bool,

    /// How many families to search at once, and the cores they share. Each
    /// search is single-threaded, so this is the whole of the parallelism.
    /// `fanout` only
    #[arg(short = 'j', long)]
    pub jobs: Option<usize>,

    /// Threads for the one search, and the cores it is pinned to. `union`
    /// only
    #[arg(short, long)]
    pub threads: Option<usize>,

    /// A recruit whose original hmmsearch hits at or below this E-value is a
    /// reject for every tool
    #[arg(short = 'e', default_value_t = REJECT_E, value_name = "F")]
    pub reject_e: f64,
}

#[derive(Parser, Debug)]
pub struct FillArgs {
    #[command(flatten)]
    pub place: Where,

    /// How to search. The two agree; they cost differently
    #[arg(long, value_enum, default_value_t = Strategy::Fanout)]
    pub strategy: Strategy,

    #[arg(long)]
    pub dry_run: bool,

    /// How many families to search at once, and the cores they share. Each
    /// search is single-threaded, so this is the whole of the parallelism.
    /// `fanout` only
    #[arg(short = 'j', long)]
    pub jobs: Option<usize>,

    /// Threads per search, and the cores each search is pinned to. `union`
    /// only, where there is one search per tool rather than one per family
    #[arg(short, long)]
    pub threads: Option<usize>,
}

#[derive(Parser, Debug)]
pub struct LearnArgs {
    #[command(flatten)]
    pub place: Where,

    /// Which layout `fill` left behind. Must match what it ran with
    #[arg(long, value_enum, default_value_t = Strategy::Fanout)]
    pub strategy: Strategy,

    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Parser, Debug)]
pub struct AllArgs {
    #[command(flatten)]
    pub place: Where,

    /// How `reject` and `fill` search. The two agree; they cost differently
    #[arg(long, value_enum, default_value_t = Strategy::Fanout)]
    pub strategy: Strategy,

    /// Threads per search, and for gather and learn the threads they read with
    #[arg(short, long)]
    pub threads: Option<usize>,

    /// How many families `reject` and `fill` search at once, and the cores
    /// they share
    #[arg(short = 'j', long)]
    pub jobs: Option<usize>,
}

/// The set this calibrates against: the same deal a benchmark searches, drawn
/// backwards. Everything it makes out of that is its own output rather than a
/// second set.
pub const SHAPE: &util::set::Shape = &util::set::shape::REVERSED;

#[derive(Parser)]
#[command(
    name = "cutoffs",
    about = "per-family score cutoffs from reversed decoys"
)]
struct Cli {
    #[command(subcommand)]
    command: Cmd,
}

/// Where this calibration reads and writes, as one label of `paths.toml` names
/// them.
///
/// `set` is what it calibrates; the rest is what it makes out of that, one
/// directory per stage.
#[derive(serde::Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct Paths {
    pub set: PathBuf,
    pub gather: PathBuf,
    pub recruit: PathBuf,
    pub reject: PathBuf,
    pub fill: PathBuf,
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
            fill: file.at(p.fill),
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
            Cmd::Fill(a) => &a.place,
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
        Cmd::Fill(args) => fill(args, paths),
        Cmd::Learn(args) => learn(args, paths),
        Cmd::All(args) => all(args, paths),
    }
}

// ----------------------------------------------------------------- recruit

fn recruit(args: RecruitArgs, paths: &Paths) -> anyhow::Result<()> {
    let threads = args.threads.context("--threads is required")?;
    let layout = Layout::new(paths)?;

    let nail_bin = nail()?;
    let mmseqs_bin = mmseqs()?;

    let query_hmm = layout.query_hmm();
    let query_db = layout.query_db();

    let stage = layout.recruit()?;
    let (results, tmp) = (stage.results(), stage.tmp());

    let mut pl = PipelineBuilder::new()
        .pool(threads)
        .step(PCmd::new("mkdir").flag("-p").path(&results));

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
                    util::search::createdb(
                        &mmseqs_bin,
                        shard,
                        &target_db,
                        &idx.to_string(),
                        threads,
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
                    .arg("-t", threads)
                    .arg("--tmp-dir", scratch.join("nail"))
                    .arg("--mmseqs-s", RECRUIT_S)
                    .arg("--seed-mode", util::search::SEED_MODE)
                    .arg("--mmseqs-max-seqs", RECRUIT_MAX_SEQS)
                    .arg("-E", util::search::EVALUE)
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
                    let cmds = util::search::Mmseqs {
                        bin: &mmseqs_bin,
                        query_db: &query_db,
                        target_db: &target_db,
                        aln_db: aln_db.clone(),
                        work: scratch.join("work"),
                        out: manifest::table_path(&results, MMSEQS, &idx.to_string()),
                        threads,
                        s: Some(RECRUIT_S.to_string()),
                        max_seqs: Some(RECRUIT_MAX_SEQS),
                        evalue: util::search::EVALUE,
                    }
                    .cmds();

                    [
                        cmds.search
                            .field(manifest::NAME, MMSEQS)
                            .field(manifest::TOOL, MMSEQS)
                            .field(manifest::SHARD, idx.to_string()),
                        cmds.convert.field(manifest::SHARD, idx.to_string()),
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

    ledger::clear(&stage.root);
    pipeline.run()?;
    ledger::record(&stage.root)
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
    //
    // fanout only. union searches the whole query set in one invocation, and
    // the split is two files per family -- 41,590 of them for Pfam -- that
    // nothing would then read
    if args.strategy == Strategy::Union {
        println!("union: leaving the query set whole");
        println!("wrote {}", layout.gather.display());
        return Ok(());
    }

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

// ------------------------------------------------------------------ reject

/// The accession inside a `db|ACC|NAME` header, if the name is one.
///
/// mmseqs reports this where nail reports the whole name. MGnify's headers are
/// a single token and the two agree; Swissprot's are not and they do not.
fn accession(name: &str) -> Option<String> {
    let mut parts = name.split('|');
    let (_db, acc) = (parts.next()?, parts.next()?);
    parts.next()?;
    (!acc.is_empty()).then(|| acc.to_string())
}

/// The family of each `<family>.fa` in a directory, sorted.
fn fa_stems(dir: &Path) -> anyhow::Result<Vec<String>> {
    let mut families: Vec<String> = std::fs::read_dir(dir)
        .with_context(|| format!("failed to read {}", dir.display()))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "fa"))
        .filter_map(|p| Some(p.file_stem()?.to_str()?.to_string()))
        .collect();
    families.sort();
    Ok(families)
}

/// The sequences each family's `<family>.fa` in a directory holds, under both
/// spellings a tool may report.
///
/// This is the map `union` needs and `fanout` gets for free from the file
/// layout: a hit only counts for the family whose file held its sequence.
fn names_by_family(dir: &Path) -> anyhow::Result<HashMap<String, HashSet<String>>> {
    let mut out: HashMap<String, HashSet<String>> = HashMap::new();

    for family in fa_stems(dir)? {
        let path = dir.join(format!("{family}.fa"));
        let names = out.entry(family).or_default();
        let mut rows = Reader::new(std::fs::File::open(&path)?, Format::Fasta);
        while rows.advance()? {
            if let Some(name) = libsail::seq::name_of(Format::Fasta, rows.record()) {
                let name = String::from_utf8_lossy(name).into_owned();
                // both spellings, because the tools do not agree on one. a
                // Swissprot record is `sp|Q7PKQ5|SQUT_ECO57`: nail reports it
                // whole and mmseqs reports the accession alone, so a map keyed
                // on either form alone silently drops one tool's every hit
                if let Some(acc) = accession(&name) {
                    names.insert(acc);
                }
                names.insert(name);
            }
        }
    }

    Ok(out)
}

/// Every record under a directory of per-family fasta once, written into one
/// file.
///
/// A sequence several families hold appears once here and is searched once;
/// `learn` is what puts each hit back with the family that recruited it.
/// Streamed rather than collected, so what is held is the set of names seen
/// and not the sequences.
fn union_fasta(from: &Path, to: &Path) -> anyhow::Result<usize> {
    let mut seen: HashSet<Vec<u8>> = HashSet::new();
    let mut out = std::io::BufWriter::with_capacity(
        1 << 20,
        std::fs::File::create(to).with_context(|| format!("failed to create {}", to.display()))?,
    );

    for family in fa_stems(from)? {
        let path = from.join(format!("{family}.fa"));
        let mut rows = Reader::new(std::fs::File::open(&path)?, Format::Fasta);
        while rows.advance()? {
            let rec = rows.record();
            let Some(name) = libsail::seq::name_of(Format::Fasta, rec) else {
                continue;
            };
            if seen.insert(name.to_vec()) {
                out.write_all(rec)?;
                if !rec.ends_with(b"\n") {
                    out.write_all(b"\n")?;
                }
            }
        }
    }

    out.flush()?;
    Ok(seen.len())
}

/// hmmsearch, as every search in this crate that is not batched through
/// `util::search::hmmer` runs it.
fn hmmsearch_cmd(bin: &Path, hmm: &Path, target: &Path, tblout: &Path, evalue: &str) -> PCmd {
    PCmd::new(bin)
        .arg("--cpu", 1)
        .arg("-E", evalue)
        .arg("-o", "/dev/null")
        .arg("--tblout", tblout)
        .path(hmm)
        .path(target)
}

fn reject(args: RejectArgs, paths: &Paths) -> anyhow::Result<()> {
    let layout = Layout::new(paths)?;
    let stage = layout.reject()?;
    let originals = layout.originals();

    if !originals.is_dir() {
        bail!(
            "no recruits in {}; run `cutoffs gather` first",
            originals.display()
        );
    }
    let families = fa_stems(&originals)?;
    if families.is_empty() {
        bail!("no recruit files in {}", originals.display());
    }

    let results = stage.results();
    if results.exists() {
        std::fs::remove_dir_all(&results)?;
    }
    let tmp = stage.tmp();
    let hmmsearch_bin = hmmsearch()?;

    let pl = match args.strategy {
        Strategy::Fanout => {
            let jobs = args.jobs.context(
                "--jobs is required: it is how many families run at once, and the cores they share",
            )?;
            println!("judging {} families, {jobs} at once...", families.len());

            let queries = layout.queries();
            PipelineBuilder::new()
                .pool(jobs)
                .step(
                    Step::serial([PCmd::new("mkdir").name("dirs").flag("-p").path(&results)])
                        .name("dirs"),
                )
                .step(
                    Step::batched(
                        jobs,
                        families.iter().map(|family| {
                            hmmsearch_cmd(
                                &hmmsearch_bin,
                                &queries.join(family).join("query.hmm"),
                                &originals.join(format!("{family}.fa")),
                                &manifest::table_path(&results, JUDGE, family),
                                util::search::EVALUE,
                            )
                            .arg("-Z", layout.seqs)
                            .field(manifest::NAME, JUDGE)
                            .field(manifest::TOOL, HMMER)
                            .field(manifest::SHARD, family)
                        }),
                    )
                    .name(JUDGE),
                )
        }
        Strategy::Union => {
            let threads = args
                .threads
                .context("--threads is required with --strategy union")?;
            if tmp.exists() {
                std::fs::remove_dir_all(&tmp)?;
            }
            std::fs::create_dir_all(&tmp)?;

            let pool = tmp.join("originals.fa");
            let n = union_fasta(&originals, &pool)?;
            println!("{n} distinct recruits in {}", pool.display());

            // hmmsearch does not scale past a couple of threads, so it gets
            // the query cut into parts and the parts run together at
            // util::search::HMMER_CPU each, rather than one invocation holding
            // the whole pool
            let split = util::search::Split::new(
                &layout.query_hmm(),
                util::split::Kind::Hmm,
                tmp.join("query"),
                util::search::jobs(threads),
            );
            let parts = split.parts();
            let scratch = tmp.join("hmmer");

            let search = parts.iter().enumerate().map(|(i, part)| {
                util::search::tag(
                    PCmd::new(&hmmsearch_bin)
                        .name(i.to_string())
                        .arg("--cpu", util::search::HMMER_CPU)
                        .arg("-Z", layout.seqs)
                        .arg("--tblout", scratch.join(format!("{i}.tbl")))
                        .arg("-E", util::search::EVALUE)
                        .path(part)
                        .path(&pool),
                    JUDGE,
                    HMMER,
                    &[],
                )
                .field(manifest::SHARD, ALL)
            });

            let cat = (0..parts.len())
                .fold(PCmd::new("cat").name("tbl"), |cmd, i| {
                    cmd.path(scratch.join(format!("{i}.tbl")))
                })
                .stdout_to(manifest::table_path(&results, JUDGE, ALL));

            PipelineBuilder::new()
                .pool(threads)
                .step(
                    Step::serial([PCmd::new("mkdir")
                        .name("dirs")
                        .flag("-p")
                        .path(&results)
                        .path(&scratch)])
                    .name("dirs"),
                )
                .step(split.step("split", &[]))
                .step(
                    Step::batched(parts.len(), search)
                        .name(JUDGE)
                        // one pool the parts share, as util::search::hmmer
                        // runs them
                        .pool(util::search::HMMER_CPU * parts.len()),
                )
                .step(
                    Step::serial([
                        util::search::tag(cat, JUDGE, HMMER, &[]).field(manifest::SHARD, ALL)
                    ])
                    .name(format!("cat.{JUDGE}")),
                )
        }
    };

    let pipeline = pl
        .stderr_dir(tmp.join("stderr"))
        .sink(Progress::new())
        .sink(PTable::new(stage.manifest()))
        .build()
        .context("failed to build the judge")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    ledger::clear(&stage.root);
    pipeline.run()?;
    ledger::record(&stage.root)?;

    // the judge's verdicts, joined to the family that recruited each
    // sequence. under union the table holds every family against every
    // original, and a hit on a sequence the family never recruited answers a
    // question nobody asked: see CLAUDE.md on rejection being per pair
    let recruits = names_by_family(&originals)?;
    let tables: Vec<PathBuf> = match args.strategy {
        Strategy::Fanout => families
            .iter()
            .map(|f| manifest::table_path(&results, JUDGE, f))
            .collect(),
        Strategy::Union => vec![manifest::table_path(&results, JUDGE, ALL)],
    };

    let mut rejects: HashMap<(String, String), (f64, f32)> = HashMap::new();
    for path in &tables {
        if !path.exists() {
            continue;
        }
        let judge = Table::<HitParser<HmmerTable>>::open(path)
            .with_context(|| format!("failed to read {}", path.display()))?;
        for h in judge.iter() {
            if h.e_value > args.reject_e {
                continue;
            }
            if !recruits
                .get(&h.query)
                .is_some_and(|names| names.contains(&h.target))
            {
                continue;
            }
            let entry = rejects
                .entry((h.query.clone(), h.target.clone()))
                .or_insert((h.e_value, h.score));
            if h.e_value < entry.0 {
                *entry = (h.e_value, h.score);
            }
        }
    }

    let mut rows: Vec<_> = rejects.into_iter().collect();
    rows.sort_by(|a, b| a.0.cmp(&b.0));
    let families_hit = rows
        .iter()
        .map(|((f, _), _)| f)
        .collect::<HashSet<_>>()
        .len();

    let mut table = toil::Table::new(toil::Schema::new(["family", "target", "evalue", "score"]));
    for ((family, target), (e, score)) in &rows {
        table.row(vec![
            family.clone(),
            target.clone(),
            format!("{e:.1e}"),
            format!("{score:.1}"),
        ]);
    }
    let out = layout.rejects_tbl();
    table
        .write(&out)
        .with_context(|| format!("failed to write {}", out.display()))?;
    println!(
        "{} rejects over {families_hit} families, written to {}",
        rows.len(),
        out.display()
    );
    Ok(())
}

// -------------------------------------------------------------------- fill

/// The rejected (family, sequence) pairs `reject` wrote.
fn read_rejects(path: &Path) -> anyhow::Result<HashSet<(String, String)>> {
    let table = toil::Table::read(path)
        .with_context(|| format!("failed to read {}", path.display()))
        .with_context(|| {
            format!(
                "no rejects at {}; run `cutoffs reject` first",
                path.display()
            )
        })?;
    let (family, target) = (
        table
            .index("family")
            .context("rejects.tbl has no family column")?,
        table
            .index("target")
            .context("rejects.tbl has no target column")?,
    );
    Ok(table
        .rows()
        .iter()
        .filter_map(|row| Some((row.get(family)?.to_string(), row.get(target)?.to_string())))
        .collect())
}

fn fill(args: FillArgs, paths: &Paths) -> anyhow::Result<()> {
    let layout = Layout::new(paths)?;
    let stage = layout.fill()?;
    let reversals = layout.reversals();

    if !reversals.is_dir() {
        bail!(
            "no recruits in {}; run `cutoffs gather` first",
            reversals.display()
        );
    }
    let rejects = read_rejects(&layout.rejects_tbl())?;

    // the decoys: every recruit the judge did not reject, laid out per family
    // the way gather laid out the recruits, and listed in decoys.tbl so that
    // learn reads what was searched rather than what is on disk
    let decoys_dir = layout.decoys();
    if decoys_dir.exists() {
        std::fs::remove_dir_all(&decoys_dir)?;
    }
    std::fs::create_dir_all(&decoys_dir)?;

    let mut listed = toil::Table::new(toil::Schema::new(["family", "target"]));
    let mut families = Vec::new();
    let (mut kept, mut dropped) = (0usize, 0usize);
    for family in fa_stems(&reversals)? {
        let mut out = std::io::BufWriter::new(std::fs::File::create(
            decoys_dir.join(format!("{family}.fa")),
        )?);
        let mut n = 0;
        let mut rows = Reader::new(
            std::fs::File::open(reversals.join(format!("{family}.fa")))?,
            Format::Fasta,
        );
        while rows.advance()? {
            let rec = rows.record();
            let Some(name) = libsail::seq::name_of(Format::Fasta, rec) else {
                continue;
            };
            let name = String::from_utf8_lossy(name).into_owned();
            if rejects.contains(&(family.clone(), name.clone())) {
                dropped += 1;
                continue;
            }
            out.write_all(rec)?;
            if !rec.ends_with(b"\n") {
                out.write_all(b"\n")?;
            }
            listed.row(vec![family.clone(), name]);
            n += 1;
        }
        out.flush()?;
        kept += n;
        if n > 0 {
            families.push(family);
        }
    }
    let decoys_tbl = layout.decoys_tbl();
    listed
        .write(&decoys_tbl)
        .with_context(|| format!("failed to write {}", decoys_tbl.display()))?;
    println!(
        "{kept} decoys over {} families, {dropped} recruits rejected",
        families.len()
    );
    if families.is_empty() {
        bail!("every recruit was rejected; nothing to fill");
    }

    let results = stage.results();
    if results.exists() {
        std::fs::remove_dir_all(&results)?;
    }
    let tmp = stage.tmp();

    let nail_bin = nail()?;
    let mmseqs_bin = mmseqs()?;
    let hmmsearch_bin = hmmsearch()?;

    let pl = match args.strategy {
        Strategy::Fanout => {
            let jobs = args.jobs.context(
                "--jobs is required: it is how many families run at once, and the cores they share",
            )?;
            println!("filling {} families, {jobs} at once...", families.len());

            let queries = layout.queries();
            let target = |family: &str| decoys_dir.join(format!("{family}.fa"));
            let hmm = |family: &str| queries.join(family).join("query.hmm");
            let table = |tool: &str, family: &str| manifest::table_path(&results, tool, family);
            let scratch = |family: &str| tmp.join(family);
            let query_db = |family: &str| scratch(family).join("queryDB");

            // one command per family per link of the chain, and one Step per
            // link. the searches here are single-query -- one family's profile
            // against its own decoys -- and every tool parallelises over
            // queries, so a thread count above 1 buys nothing and the
            // parallelism has to come from running many families at once.
            // `batched` is that: `jobs` commands in flight, each of them
            // single-threaded.
            //
            // which is why the steps are links rather than families. a
            // family's commands have to run in order, and a Step holds Cmds
            // rather than Steps, so a batch of ordered chains is not a thing
            // michi can be asked for. the transpose is: every family's link k,
            // then every family's link k+1. the ordering each family needs
            // still holds, since a step finishes before the next one starts.
            let per_family = |f: &dyn Fn(&str) -> PCmd| -> Vec<PCmd> {
                families.iter().map(|family| f(family)).collect()
            };

            let mmseqs_cmds = |family: &str| {
                let d = scratch(family);
                util::search::Mmseqs {
                    bin: &mmseqs_bin,
                    query_db: &query_db(family),
                    target_db: &d.join("targetDB/targetDB"),
                    aln_db: d.join("alnDB/alnDB"),
                    work: d.join("work"),
                    out: table(MMSEQS, family),
                    threads: 1,
                    s: Some(FILL_S.to_string()),
                    max_seqs: Some(FILL_MAX_SEQS),
                    evalue: FILL_E,
                }
                .cmds()
            };

            // every batch below runs in one pool of `jobs` cores, which its
            // commands share rather than each leasing one
            PipelineBuilder::new()
                .pool(jobs)
                .step(
                    Step::batched(
                        jobs,
                        per_family(&|family| {
                            let d = scratch(family);
                            PCmd::new("mkdir")
                                .name("dirs")
                                .flag("-p")
                                .path(&results)
                                .path(d.join("targetDB"))
                                .path(d.join("alnDB"))
                                .field(manifest::STAGE, DIRS)
                        }),
                    )
                    .name("dirs"),
                )
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
                )
                .step(
                    Step::batched(
                        jobs,
                        per_family(&|family| {
                            PCmd::new(&nail_bin)
                                .sub("search")
                                .arg("--mmseqs-path", &mmseqs_bin)
                                .arg("-t", 1)
                                .arg("--tmp-dir", scratch(family).join("nail"))
                                .arg("--mmseqs-s", FILL_S)
                                .arg("--seed-mode", util::search::SEED_MODE)
                                .arg("--mmseqs-max-seqs", FILL_MAX_SEQS)
                                .arg("-E", FILL_E)
                                .flag("--allow-overwrite")
                                .arg("--tbl-out", table(NAIL, family))
                                .path(hmm(family))
                                .path(target(family))
                                .field(manifest::NAME, NAIL)
                                .field(manifest::TOOL, NAIL)
                                .field(manifest::SHARD, family)
                        }),
                    )
                    .name(NAIL),
                )
                .step(
                    Step::batched(
                        jobs,
                        per_family(&|family| {
                            util::search::createdb(
                                &mmseqs_bin,
                                &target(family),
                                &scratch(family).join("targetDB/targetDB"),
                                family,
                                1,
                            )
                        }),
                    )
                    .name("createdb"),
                )
                .step(
                    Step::batched(
                        jobs,
                        per_family(&|family| {
                            mmseqs_cmds(family)
                                .search
                                .field(manifest::NAME, MMSEQS)
                                .field(manifest::TOOL, MMSEQS)
                                .field(manifest::SHARD, family)
                        }),
                    )
                    .name(MMSEQS),
                )
                .step(
                    Step::batched(
                        jobs,
                        per_family(&|family| {
                            mmseqs_cmds(family).convert.field(manifest::SHARD, family)
                        }),
                    )
                    .name("convert"),
                )
                .step(
                    Step::batched(
                        jobs,
                        // one family per command, batched over families, so
                        // there is no split to cut and nothing to concatenate
                        // the way util::search::hmmer does
                        per_family(&|family| {
                            hmmsearch_cmd(
                                &hmmsearch_bin,
                                &hmm(family),
                                &target(family),
                                &table(HMMER, family),
                                FILL_E,
                            )
                            .arg("--domtblout", manifest::dom_path(&results, HMMER, family))
                            .field(manifest::NAME, HMMER)
                            .field(manifest::TOOL, HMMER)
                            .field(manifest::SHARD, family)
                        }),
                    )
                    .name(HMMER),
                )
                .step(
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
                )
        }
        Strategy::Union => {
            let threads = args
                .threads
                .context("--threads is required with --strategy union")?;
            if tmp.exists() {
                std::fs::remove_dir_all(&tmp)?;
            }
            std::fs::create_dir_all(&tmp)?;

            let pool = tmp.join("decoys.fa");
            let n = union_fasta(&decoys_dir, &pool)?;
            println!("{n} distinct decoys in {}", pool.display());

            let query_hmm = layout.query_hmm();
            let query_db = layout.query_db();
            let table = |tool: &str| manifest::table_path(&results, tool, ALL);
            let scratch = tmp.join("search");

            let split = util::search::Split::new(
                &query_hmm,
                util::split::Kind::Hmm,
                tmp.join("query"),
                util::search::jobs(threads),
            );

            let cmds = util::search::Mmseqs {
                bin: &mmseqs_bin,
                query_db: &query_db,
                target_db: &scratch.join("targetDB/targetDB"),
                aln_db: scratch.join("alnDB/alnDB"),
                work: scratch.join("work"),
                out: table(MMSEQS),
                threads,
                s: Some(FILL_S.to_string()),
                max_seqs: Some(FILL_MAX_SEQS),
                evalue: FILL_E,
            }
            .cmds();

            let hmmer = util::search::hmmer(
                &hmmsearch_bin,
                &split,
                &util::search::Dirs::new(&stage.root, &scratch),
                HMMER,
                HMMER,
                ALL,
                &pool,
                FILL_E,
                true,
                &[],
            );

            PipelineBuilder::new()
                .pool(threads)
                .step(
                    Step::serial([PCmd::new("mkdir")
                        .name("dirs")
                        .flag("-p")
                        .path(&results)
                        .path(scratch.join("targetDB"))
                        .path(scratch.join("alnDB"))
                        .path(scratch.join("hmmer"))])
                    .name("dirs"),
                )
                .step(split.step("split", &[]))
                .step(
                    Step::serial([PCmd::new(&nail_bin)
                        .sub("search")
                        .arg("--mmseqs-path", &mmseqs_bin)
                        .arg("-t", threads)
                        .arg("--tmp-dir", scratch.join("nail"))
                        .arg("--mmseqs-s", FILL_S)
                        .arg("--seed-mode", util::search::SEED_MODE)
                        .arg("--mmseqs-max-seqs", FILL_MAX_SEQS)
                        .arg("-E", FILL_E)
                        .flag("--allow-overwrite")
                        .arg("--tbl-out", table(NAIL))
                        .path(&query_hmm)
                        .path(&pool)
                        .field(manifest::NAME, NAIL)
                        .field(manifest::TOOL, NAIL)
                        .field(manifest::SHARD, ALL)])
                    .name(NAIL)
                    .cores(threads),
                )
                .step(
                    Step::serial([util::search::createdb(
                        &mmseqs_bin,
                        &pool,
                        &scratch.join("targetDB/targetDB"),
                        ALL,
                        threads,
                    )])
                    .name("createdb"),
                )
                .step(
                    Step::serial([cmds
                        .search
                        .field(manifest::NAME, MMSEQS)
                        .field(manifest::TOOL, MMSEQS)
                        .field(manifest::SHARD, ALL)])
                    .name(MMSEQS)
                    .cores(threads),
                )
                .step(Step::serial([cmds.convert.field(manifest::SHARD, ALL)]).name("convert"))
                .step(hmmer.search.name(HMMER))
                .step(hmmer.cat.name(format!("cat.{HMMER}")))
        }
    };

    let pipeline = pl
        .stderr_dir(tmp.join("stderr"))
        .sink(Progress::new())
        .sink(PTable::new(stage.manifest()))
        .build()
        .context("failed to build the fill")?;

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

/// The decoys each family was filled against, under both spellings a tool
/// may report.
fn decoys_by_family(path: &Path) -> anyhow::Result<HashMap<String, HashSet<String>>> {
    let table = toil::Table::read(path)
        .with_context(|| format!("failed to read {}", path.display()))
        .with_context(|| format!("no decoys at {}; run `cutoffs fill` first", path.display()))?;
    let (family, target) = (
        table
            .index("family")
            .context("decoys.tbl has no family column")?,
        table
            .index("target")
            .context("decoys.tbl has no target column")?,
    );
    let mut out: HashMap<String, HashSet<String>> = HashMap::new();
    for row in table.rows() {
        let (Some(f), Some(t)) = (row.get(family), row.get(target)) else {
            continue;
        };
        let names = out.entry(f.to_string()).or_default();
        if let Some(acc) = accession(t) {
            names.insert(acc);
        }
        names.insert(t.to_string());
    }
    Ok(out)
}

/// One tool's table read into the top `N_SCORES` decoy scores and the decoy
/// count of every family it holds.
///
/// A hit counts for a family only if that family was filled against the
/// sequence: a family's null holds only sequences vetted for that family. See
/// CLAUDE.md. Under fanout the table is one family's and the check is moot;
/// under union it is the whole point.
fn decoy_scores<T>(
    path: &Path,
    decoys: &HashMap<String, HashSet<String>>,
) -> anyhow::Result<HashMap<String, (Vec<f32>, usize)>>
where
    T: HitColumns,
{
    let table = Table::<HitParser<T>>::open(path)
        .with_context(|| format!("failed to read {}", path.display()))?;

    // one entry per pair, the last row for it winning, so a pair reported
    // twice counts as one decoy
    let mut by_pair: HashMap<(&str, &str), &Hit> = HashMap::new();
    for hit in table.iter() {
        let pair = (hit.query.as_str(), hit.target.as_str());
        if !decoys
            .get(pair.0)
            .is_some_and(|names| names.contains(pair.1))
        {
            continue;
        }
        by_pair.insert(pair, hit);
    }

    let mut by_family: HashMap<&str, Vec<&Hit>> = HashMap::new();
    for ((family, _), hit) in by_pair {
        by_family.entry(family).or_default().push(hit);
    }

    let mut out = HashMap::new();
    for (family, mut hits) in by_family {
        hits.sort_by(|a, b| {
            a.e_value
                .partial_cmp(&b.e_value)
                .expect("NaN in decoy e-values")
        });
        let scores: Vec<f32> = hits
            .iter()
            .map(|h| h.score)
            // a family with fewer than N_SCORES decoys pads with zero, which
            // the readers of cutoffs.tbl take as "no usable cutoff"
            .chain(std::iter::repeat(0.0))
            .take(N_SCORES)
            .collect();
        out.insert(family.to_string(), (scores, hits.len()));
    }
    Ok(out)
}

fn learn(args: LearnArgs, paths: &Paths) -> anyhow::Result<()> {
    let layout = Layout::new(paths)?;
    let stage = layout.fill()?;
    let results = stage.results();

    // what the fill actually did, rather than what is on disk: a run that
    // died leaves a table behind, and a half-written one reads as a family
    // with fewer decoys than it has
    let searched = manifest::Manifest::read(&stage.manifest())?;
    let failed = searched.failed().count();
    if failed > 0 {
        bail!(
            "{failed} searches in {} did not finish; re-run `cutoffs fill`",
            stage.manifest().display()
        );
    }

    let decoys = decoys_by_family(&layout.decoys_tbl())?;
    let mut families: Vec<String> = decoys.keys().cloned().collect();
    families.sort_unstable();
    if families.is_empty() {
        bail!("no decoys in {}", layout.decoys_tbl().display());
    }

    let empty = || (vec![0.0; N_SCORES], 0);

    let rows: Vec<Vec<String>> = match args.strategy {
        Strategy::Union => {
            // read once for every family, rather than once per family
            let table = |tool: &str| manifest::table_path(&results, tool, ALL);
            let nail = decoy_scores::<NailTable>(&table(NAIL), &decoys)?;
            let mmseqs = decoy_scores::<BlastTable>(&table(MMSEQS), &decoys)?;
            let hmmer = decoy_scores::<HmmerTable>(&table(HMMER), &decoys)?;
            families
                .iter()
                .map(|family| {
                    let mut row = vec![family.clone()];
                    for scores in [&nail, &mmseqs, &hmmer] {
                        row.extend(cells(&scores.get(family).cloned().unwrap_or_else(empty)));
                    }
                    row
                })
                .collect()
        }
        Strategy::Fanout => {
            let skipped = AtomicUsize::new(0);
            // collected rather than written as they finish, so the file is in
            // family order however the pool interleaves
            let pool = pool(args.threads)?;
            let rows: Vec<Option<Vec<String>>> = pool.install(|| {
                families
                    .par_iter()
                    .map(|family| -> anyhow::Result<Option<Vec<String>>> {
                        let table = |tool: &str| manifest::table_path(&results, tool, family);
                        let (nail, mmseqs, hmmer) = (table(NAIL), table(MMSEQS), table(HMMER));
                        if !nail.exists() || !mmseqs.exists() {
                            // a family without both tables tells us nothing
                            // comparative
                            skipped.fetch_add(1, Ordering::Relaxed);
                            return Ok(None);
                        }
                        let one = |scores: HashMap<String, (Vec<f32>, usize)>| {
                            scores.get(family.as_str()).cloned().unwrap_or_else(empty)
                        };
                        let mut row = vec![family.clone()];
                        row.extend(cells(&one(decoy_scores::<NailTable>(&nail, &decoys)?)));
                        row.extend(cells(&one(decoy_scores::<BlastTable>(&mmseqs, &decoys)?)));
                        row.extend(cells(&match hmmer.exists() {
                            true => one(decoy_scores::<HmmerTable>(&hmmer, &decoys)?),
                            false => empty(),
                        }));
                        Ok(Some(row))
                    })
                    .collect::<anyhow::Result<Vec<_>>>()
            })?;
            let skipped = skipped.load(Ordering::Relaxed);
            if skipped > 0 {
                eprintln!("skipped {skipped} families missing a nail or mmseqs table");
            }
            rows.into_iter().flatten().collect()
        }
    };

    let mut headers = vec!["family".to_string()];
    for tool in TOOLS {
        headers.extend((1..=N_SCORES).map(|i| format!("{tool}_{i}")));
        headers.push(format!("{tool}_n"));
    }

    let mut table = toil::Table::new(toil::Schema::new(headers));
    for row in rows {
        table.row(row);
    }

    let out_path = layout.cutoffs_tbl()?;
    table
        .write(&out_path)
        .with_context(|| format!("failed to write {}", out_path.display()))?;
    println!("wrote {}", out_path.display());
    Ok(())
}

/// One tool's cells of a `cutoffs.tbl` row: the scores, then the count.
fn cells(scores: &(Vec<f32>, usize)) -> Vec<String> {
    scores
        .0
        .iter()
        .map(|s| format!("{s:.1}"))
        .chain(std::iter::once(scores.1.to_string()))
        .collect()
}

// --------------------------------------------------------------------- all

fn all(args: AllArgs, paths: &Paths) -> anyhow::Result<()> {
    let threads = args.threads.context("--threads is required")?;

    recruit(
        RecruitArgs {
            place: args.place.clone(),
            threads: Some(threads),
            dry_run: false,
        },
        paths,
    )?;

    gather(
        GatherArgs {
            place: args.place.clone(),
            strategy: args.strategy,
            threads,
        },
        paths,
    )?;

    reject(
        RejectArgs {
            place: args.place.clone(),
            strategy: args.strategy,
            threads: Some(threads),
            dry_run: false,
            jobs: args.jobs,
            reject_e: REJECT_E,
        },
        paths,
    )?;

    fill(
        FillArgs {
            place: args.place.clone(),
            strategy: args.strategy,
            threads: Some(threads),
            dry_run: false,
            jobs: args.jobs,
        },
        paths,
    )?;

    learn(
        LearnArgs {
            place: args.place.clone(),
            strategy: args.strategy,
            threads,
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

#[cfg(test)]
mod tests {
    use super::*;

    /// nail reports a Swissprot target whole and mmseqs reports the accession
    /// alone, so the recruit map has to answer to both. Keyed on one form
    /// only, the other tool's every hit is dropped and its cutoffs come out
    /// zero without anything failing.
    #[test]
    fn an_accession_is_read_out_of_a_piped_name() {
        assert_eq!(accession("sp|Q7PKQ5|SQUT_ECO57").as_deref(), Some("Q7PKQ5"));
        assert_eq!(accession("tr|A0A1B2|SOME_NAME").as_deref(), Some("A0A1B2"));
    }

    /// MGnify's headers are one token, which is why this went unnoticed: there
    /// the two tools agree and there is nothing to strip.
    #[test]
    fn a_plain_name_has_no_accession_inside_it() {
        assert_eq!(accession("MGYP001482868479"), None);
        assert_eq!(accession("sp|Q7PKQ5"), None);
        assert_eq!(accession(""), None);
    }
}
