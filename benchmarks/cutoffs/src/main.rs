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
//!   1. `recruit` — a sweep of nail, mmseqs and hmmer over the initial
//!      reversals, which finds the small subset that scores at all. What it
//!      finds are *recruits*.
//!   2. `gather` — pulls each recruit out of the shard it came from. A record
//!      read there is already the reversal, and reversing it again gives the
//!      *original*, so both forms come out of one pass. `recruits.tbl` lists
//!      every recruit and the tools that recruited it.
//!   3. `reject` — one hmmsearch of each family against its recruits'
//!      originals, at the `-Z` of one shard, or one pass over the tables of a
//!      forward search of the same sequences. A recruit whose original the
//!      family hits is a *reject*, a reversed homolog rather than a piece of
//!      noise; the rest are *decoys*.
//!   4. `fill` — nail and mmseqs over the decoys each did not recruit, with
//!      the prefilter effectively off and `-E 1e9`. It rescores nothing: a
//!      tool's score for a pair does not move with these settings, so a decoy
//!      a tool recruited keeps the score recruit gave it. hmmer has no gaps,
//!      since its recruit search ran over every pair in each shard.
//!   5. `learn` — the top scores per family per tool, from recruit and fill,
//!      into `cutoffs.tbl`.
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
//! straight into steps. `reject` and `fill` are the opposite shape: one
//! search per family, and every tool here parallelises over queries alone --
//! so a thread count above one buys nothing and the parallelism has to be
//! many families at once, which is a batched step of single-threaded
//! commands. Each step is one tool over the families it has work in;
//! mmseqs' five-process chain runs as one shell command, so michi launches
//! it once.
//!
//! Every launch forks this process, so what it holds in memory while a step
//! runs is what each of those 20,000 launches pays for. The stages write
//! their tables as streams and make and remove directories in this thread,
//! and nothing short-lived is a command of its own.
//!
//! `gather` and `learn` run no tools, so they stay a plain rayon pool.

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

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
use libsail::tbl::{HitColumns, HitParser, Table};
use michi::{Cmd as PCmd, PipelineBuilder, Progress, Step, Table as PTable};
use util::search::{Dirs, Split};
use util::split::Kind;
use util::tools::{hmmsearch, mmseqs, nail};
use util::{ledger, manifest};

use util::set::Set;

use util::cut;

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

// wide open, so every decoy a tool can align gets a row: the point of the
// fill is a score from every tool for every decoy, including the ones that
// tool never recruited, and a reporting threshold would hide the weak ones.
// one E-value for the three tools rather than a -Z for the two that take one
const FILL_S: &str = "12.0";
const FILL_MAX_SEQS: usize = 1_000_000_000;
const FILL_E: &str = "1e9";

// static rather than the shared prog mode: prog grows a query's seed list in
// rounds and stops it by the hit rate of the last round, so on a target of
// two sequences and on a pool of two thousand it reaches different pairs.
// static aligns the whole prefilter list, which on these targets is every
// decoy, and that is what the fill is for
const FILL_SEED_MODE: &str = "static";

// nail hands its mmseqs align step `-e S * Z`, with Z the target count, and
// mmseqs scales its E-value by residues rather than by sequences, so the
// effective filter moves with the mean length of the target file and two
// families' decoy files get different ones. a P-value no hit exceeds turns
// the step off, the way FILL_E does for the reporting threshold
const FILL_SEED_P: &str = "1e9";

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

/// The column width `rejects.tbl` and `decoys.tbl` are streamed at for a
/// sequence name; a longer one overruns its column and still reads.
const NAME_WIDTH: usize = 32;

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
    // small file of one family's recruits. a deal's shards differ by one
    // sequence when the count does not divide, so this is the larger, and
    // `new` refuses a set whose shards differ by more
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

        let (mut seqs, mut least) = (0u64, u64::MAX);
        for unit in set.units() {
            let n = unit.number("seqs")?;
            seqs = seqs.max(n);
            least = least.min(n);
        }
        if seqs - least > 1 {
            bail!(
                "the set at {} has shards of different sizes ({least} and {seqs} sequences), so there is no one -Z for the judge",
                paths.set.display()
            );
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

    /// Every recruit and the tools that recruited it, written by `gather`.
    fn recruits_tbl(&self) -> PathBuf {
        self.gather.join("recruits.tbl")
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

    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Parser, Debug)]
pub struct RejectArgs {
    #[command(flatten)]
    pub place: Where,

    #[arg(long)]
    pub dry_run: bool,

    /// How many families to search at once, and the cores they share. Each
    /// search is single-threaded, so this is the whole of the parallelism
    #[arg(short = 'j', long)]
    pub jobs: usize,

    /// A recruit whose original hmmsearch hits at or below this E-value is a
    /// reject for every tool
    #[arg(short = 'e', default_value_t = REJECT_E, value_name = "F")]
    pub reject_e: f64,

    /// Judge by reading hmmsearch's tables of a forward search of the same
    /// sequences, `hmmer.<unit>.tbl` per unit of --reject-set, rather than
    /// searching
    #[arg(long, value_name = "dir", requires = "reject_set")]
    pub reject_with: Option<PathBuf>,

    /// The set those tables searched, for each unit's sequence count
    #[arg(long, value_name = "dir", requires = "reject_with")]
    pub reject_set: Option<PathBuf>,
}

#[derive(Parser, Debug)]
pub struct FillArgs {
    #[command(flatten)]
    pub place: Where,

    #[arg(long)]
    pub dry_run: bool,

    /// How many families to search at once, and the cores they share. Each
    /// search is single-threaded, so this is the whole of the parallelism
    #[arg(short = 'j', long)]
    pub jobs: usize,
}

#[derive(Parser, Debug)]
pub struct LearnArgs {
    #[command(flatten)]
    pub place: Where,

    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Parser, Debug)]
pub struct AllArgs {
    #[command(flatten)]
    pub place: Where,

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
    let hmmsearch_bin = hmmsearch()?;

    let query_hmm = layout.query_hmm();
    let query_db = layout.query_db();

    let stage = layout.recruit()?;
    let (results, tmp) = (stage.results(), stage.tmp());
    let dirs = Dirs::new(&stage.root, &tmp);

    // one split for every shard, outside the per-shard scratch that each
    // shard's last step removes
    let split = Split::new(
        &query_hmm,
        Kind::Hmm,
        tmp.join("hmmer-query"),
        util::search::jobs(threads),
    );

    let mut pl = PipelineBuilder::new()
        .pool(threads)
        .step(PCmd::new("mkdir").flag("-p").path(&results))
        .step(split.step("split", &[]));

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
            );

        let hmmer = util::search::hmmer(
            &hmmsearch_bin,
            &split,
            &dirs,
            HMMER,
            HMMER,
            &idx.to_string(),
            shard,
            util::search::EVALUE,
            false,
            &[],
        );
        pl = pl
            .step(hmmer.search)
            .step(hmmer.cat)
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

    // per shard, which families each recruited sequence hit and which tools
    // reported it. Only names travel here; the sequences themselves are read
    // in the second pass.
    let wanted: Vec<(String, Recruits)> = pool.install(|| {
        shard_list
            .par_iter()
            .map(|shard| -> anyhow::Result<(String, Recruits)> {
                let mut map = Recruits::new();
                let path = |tool: &str| manifest::table_path(&recruit_results, tool, shard);
                let read = |path: PathBuf| format!("failed to read {}", path.display());

                let tbl = Table::<HitParser<NailTable>>::open(path(NAIL))
                    .with_context(|| read(path(NAIL)))?;
                collect(&tbl, bit(NAIL), &mut map);
                let tbl = Table::<HitParser<BlastTable>>::open(path(MMSEQS))
                    .with_context(|| read(path(MMSEQS)))?;
                collect(&tbl, bit(MMSEQS), &mut map);
                let tbl = Table::<HitParser<HmmerTable>>::open(path(HMMER))
                    .with_context(|| read(path(HMMER)))?;
                collect(&tbl, bit(HMMER), &mut map);

                Ok((shard.clone(), map))
            })
            .collect::<anyhow::Result<Vec<_>>>()
    })?;

    let families: HashSet<String> = wanted.iter().flat_map(|(_, m)| m.keys().cloned()).collect();

    if families.is_empty() {
        bail!("no family recruited any decoys; is the recruit stage's E-value too strict?");
    }

    println!("{} families recruited decoys", families.len());

    // every recruit and the tools that reported it, sorted by family and
    // then name, for reject and fill to read rather than recruit's tables
    let mut sorted: Vec<&String> = families.iter().collect();
    sorted.sort_unstable();
    std::fs::create_dir_all(&layout.gather)?;
    let mut table = recruits_stream(&layout.recruits_tbl(), &families)?;
    let mut recruits = 0usize;
    for family in sorted {
        let mut rows: Vec<(&str, u8)> = wanted
            .iter()
            .filter_map(|(_, map)| map.get(family))
            .flat_map(|names| names.iter().map(|(n, b)| (n.as_str(), *b)))
            .collect();
        rows.sort_unstable();
        for (name, bits) in rows {
            let mut row = vec![family.clone(), name.to_string()];
            row.extend(TOOLS.map(|tool| u8::from(bits & bit(tool) != 0).to_string()));
            table.row(row)?;
            recruits += 1;
        }
    }
    table.header()?;
    drop(table);
    println!(
        "{recruits} recruits, listed in {}",
        layout.recruits_tbl().display()
    );

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
                    for name in names.keys() {
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

/// Each recruited sequence of one shard, by family, with the tools that
/// reported it as bits.
type Recruits = HashMap<String, HashMap<String, u8>>;

/// The bit that marks a pair as reported by `tool` in recruit.
fn bit(tool: &str) -> u8 {
    let at = TOOLS
        .iter()
        .position(|t| *t == tool)
        .expect("a tool the calibration scores");
    1 << at
}

/// Fold one tool's hit table into a shard's recruits.
fn collect<C: HitColumns>(tbl: &Table<HitParser<C>>, tool: u8, map: &mut Recruits) {
    for hit in tbl.iter() {
        *map.entry(hit.query.clone())
            .or_default()
            .entry(hit.target.clone())
            .or_default() |= tool;
    }
}

/// `recruits.tbl` opened for writing: a family, a sequence, and a column
/// per tool holding 1 where that tool recruited the pair.
fn recruits_stream(
    path: &Path,
    families: &HashSet<String>,
) -> anyhow::Result<toil::Stream<std::io::BufWriter<std::fs::File>>> {
    let mut columns = vec!["family", "target"];
    columns.extend(TOOLS);
    let schema = toil::Schema::new(columns);
    let mut widths = schema.widths();
    let mut wide = vec![
        "x".repeat(families.iter().map(String::len).max().unwrap_or(0)),
        "x".repeat(NAME_WIDTH),
    ];
    wide.extend(TOOLS.map(|_| "0".to_string()));
    widths.widen(&schema.row(wide));
    Ok(toil::Stream::new(
        schema,
        widths,
        std::io::BufWriter::with_capacity(
            1 << 20,
            std::fs::File::create(path)
                .with_context(|| format!("failed to create {}", path.display()))?,
        ),
    ))
}

// ------------------------------------------------------------------ reject

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

/// hmmsearch on one thread, as the judge and the fill run it.
fn hmmsearch_cmd(bin: &Path, hmm: &Path, target: &Path, tblout: &Path, evalue: &str) -> PCmd {
    PCmd::new(bin)
        .arg("--cpu", 1)
        .arg("-E", evalue)
        .arg("-o", "/dev/null")
        .arg("--tblout", tblout)
        .path(hmm)
        .path(target)
}

/// One step of `jobs` single-threaded commands, one per family.
//
// the searches here are single-query -- one family's profile against its
// own recruits -- and every tool parallelises over queries, so a thread count
// above 1 buys nothing and the parallelism has to come from running many
// families at once. `batched` is that: `jobs` commands in flight, each of
// them single-threaded, sharing one pool of `jobs` cores
fn per_family(
    name: &str,
    jobs: usize,
    families: &[String],
    cmd: impl Fn(&str) -> PCmd,
) -> Step<'static> {
    Step::batched(jobs, families.iter().map(|f| cmd(f))).name(name)
}

fn reject(args: RejectArgs, paths: &Paths) -> anyhow::Result<()> {
    let layout = Layout::new(paths)?;
    let stage = layout.reject()?;
    let originals = layout.originals();

    if let (Some(tables), Some(set)) = (&args.reject_with, &args.reject_set) {
        return reject_with(
            &layout,
            &stage,
            &paths.set,
            tables,
            set,
            args.reject_e,
            args.jobs,
        );
    }

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
    std::fs::create_dir_all(&results)?;
    let tmp = stage.tmp();
    let hmmsearch_bin = hmmsearch()?;
    let queries = layout.queries();

    println!(
        "judging {} families, {} at once...",
        families.len(),
        args.jobs
    );

    let pipeline = PipelineBuilder::new()
        .pool(args.jobs)
        .step(per_family(JUDGE, args.jobs, &families, |family| {
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
        }))
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

    // the judge's verdicts. each table is one family's profile against that
    // family's own recruits, so every hit in it is a question that family
    // was asked: see CLAUDE.md on rejection being per pair
    let out = layout.rejects_tbl();
    let mut table = rejects_stream(&out, families.iter().map(String::len).max())?;

    let (mut rejects, mut families_hit) = (0usize, 0usize);
    for family in &families {
        let path = manifest::table_path(&results, JUDGE, family);
        if !path.exists() {
            continue;
        }
        let judge = Table::<HitParser<HmmerTable>>::open(&path)
            .with_context(|| format!("failed to read {}", path.display()))?;
        // one row per pair, at its best E-value
        let mut best: HashMap<&str, (f64, f32)> = HashMap::new();
        for h in judge.iter() {
            if h.e_value > args.reject_e {
                continue;
            }
            let entry = best
                .entry(h.target.as_str())
                .or_insert((h.e_value, h.score));
            if h.e_value < entry.0 {
                *entry = (h.e_value, h.score);
            }
        }
        if best.is_empty() {
            continue;
        }
        families_hit += 1;
        let mut rows: Vec<_> = best.into_iter().collect();
        rows.sort_by(|a, b| a.0.cmp(b.0));
        for (target, (e, score)) in rows {
            table.row([
                family.clone(),
                target.to_string(),
                format!("{e:.1e}"),
                format!("{score:.1}"),
            ])?;
            rejects += 1;
        }
    }
    table.header()?;
    println!(
        "{rejects} rejects over {families_hit} families, written to {}",
        out.display()
    );
    Ok(())
}

/// `rejects.tbl` opened for writing, its family column `family_width` wide.
fn rejects_stream(
    path: &Path,
    family_width: Option<usize>,
) -> anyhow::Result<toil::Stream<std::io::BufWriter<std::fs::File>>> {
    let schema = toil::Schema::new(["family", "target", "evalue", "score"]);
    let mut widths = schema.widths();
    widths.widen(&schema.row([
        "x".repeat(family_width.unwrap_or(0)),
        "x".repeat(NAME_WIDTH),
        "0.0e-000".to_string(),
        "0000.0".to_string(),
    ]));
    Ok(toil::Stream::new(
        schema,
        widths,
        std::io::BufWriter::new(
            std::fs::File::create(path)
                .with_context(|| format!("failed to create {}", path.display()))?,
        ),
    ))
}

/// A (family, sequence) pair as one number, for a set too large to hold as
/// strings.
//
// two 64-bit hashes of the pair under different leading tags. two pairs
// colliding in both is ~1e-20 at 600 million pairs, and a collision could
// only turn a forward hit on a sequence that is not a recruit into a reject
fn pair_key(family: &[u8], target: &[u8]) -> u128 {
    use std::hash::{DefaultHasher, Hash, Hasher};
    let half = |tag: u8| {
        let mut h = DefaultHasher::new();
        (tag, family, target).hash(&mut h);
        h.finish()
    };
    (u128::from(half(0)) << 64) | u128::from(half(1))
}

/// Judge every recruit in one pass over a forward search's hmmer tables.
///
/// The forward set holds the originals under the names the reversals carry,
/// dealt into shards of its own, so a recruit is found by name whichever
/// shard holds it.
fn reject_with(
    layout: &Layout,
    stage: &Stage,
    reversals: &Path,
    tables: &Path,
    set: &Path,
    reject_e: f64,
    threads: usize,
) -> anyhow::Result<()> {
    let forward = Set::load_as(set, &util::set::shape::FIXED)?;
    let reversed = Set::load_as(reversals, SHAPE)?;
    let total = |set: &Set| -> anyhow::Result<u64> { set.units().map(|u| u.number("seqs")).sum() };
    let (f, r) = (total(&forward)?, total(&reversed)?);
    if f != r {
        bail!(
            "{} holds {f} sequences and {} holds {r}; the forward search is of a different set",
            set.display(),
            reversals.display()
        );
    }

    // each table's E-values rescaled from its own shard's size to the
    // judge's -Z, so a reject means what the judge's search would make it.
    // a shard within one sequence of it is the same -Z, as Layout takes a
    // deal's shards to be: hmmer prints two figures, and rescaling those by
    // a part in a million moves a pair printed at the threshold across it
    let units = forward
        .units()
        .map(|u| {
            let path = tables.join(format!("{HMMER}.{}.tbl", u.name()));
            if !path.exists() {
                bail!("no table for unit {:?} at {}", u.name(), path.display());
            }
            let n = u.number("seqs")?;
            let scale = match n.abs_diff(layout.seqs) <= 1 {
                true => 1.0,
                false => layout.seqs as f64 / n as f64,
            };
            Ok((path, scale))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;

    let recruits_tbl = layout.recruits_tbl();
    println!("reading {}...", recruits_tbl.display());
    let mut keys: Vec<u128> = Vec::new();
    let mut family_width = 0;
    let mut rows = toil::Reader::new(std::io::BufReader::with_capacity(
        1 << 20,
        std::fs::File::open(&recruits_tbl).with_context(|| {
            format!(
                "no recruits at {}; run `cutoffs gather` first",
                recruits_tbl.display()
            )
        })?,
    ))?;
    while let Some((_, entry)) = rows.next_entry()? {
        let toil::Entry::Row(cells) = entry else {
            continue;
        };
        let (Some(family), Some(target)) = (cells.field(0), cells.field(1)) else {
            continue;
        };
        family_width = family_width.max(family.len());
        keys.push(pair_key(family, target));
    }
    keys.sort_unstable();
    keys.dedup();

    let results = stage.results();
    if results.exists() {
        std::fs::remove_dir_all(&results)?;
    }
    std::fs::remove_file(stage.manifest()).ok();
    ledger::clear(&stage.root);
    std::fs::create_dir_all(&stage.root)?;

    println!(
        "judging {} recruits against {} tables, {threads} at once...",
        keys.len(),
        units.len()
    );
    let pool = pool(threads)?;
    let mut found: Vec<(String, String, f64, f32)> = pool
        .install(|| {
            units
                .par_iter()
                .map(|(path, scale)| -> anyhow::Result<Vec<_>> {
                    let tbl = Table::<HitParser<HmmerTable>>::open(path)
                        .with_context(|| format!("failed to read {}", path.display()))?;
                    let mut hits = Vec::new();
                    for h in tbl.iter() {
                        let e = h.e_value * scale;
                        if e > reject_e {
                            continue;
                        }
                        let key = pair_key(h.query.as_bytes(), h.target.as_bytes());
                        if keys.binary_search(&key).is_ok() {
                            hits.push((h.query.clone(), h.target.clone(), e, h.score));
                        }
                    }
                    Ok(hits)
                })
                .collect::<anyhow::Result<Vec<_>>>()
        })?
        .into_iter()
        .flatten()
        .collect();
    drop(keys);

    // one row per pair at its best E-value, as the judge writes them
    found.sort_by(|a, b| {
        (&a.0, &a.1)
            .cmp(&(&b.0, &b.1))
            .then(a.2.partial_cmp(&b.2).expect("NaN in an E-value"))
    });
    found.dedup_by(|b, a| a.0 == b.0 && a.1 == b.1);

    let out = layout.rejects_tbl();
    let mut table = rejects_stream(&out, Some(family_width))?;
    table.meta("judged-by", [tables.display().to_string()])?;
    let mut families = HashSet::new();
    for (family, target, e, score) in &found {
        families.insert(family.as_str());
        table.row([
            family.clone(),
            target.clone(),
            format!("{e:.1e}"),
            format!("{score:.1}"),
        ])?;
    }
    table.header()?;
    println!(
        "{} rejects over {} families, written to {}",
        found.len(),
        families.len(),
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
    let all_families = fa_stems(&reversals)?;

    let results = stage.results();
    if results.exists() {
        std::fs::remove_dir_all(&results)?;
    }
    std::fs::create_dir_all(&results)?;
    let tmp = stage.tmp();
    if tmp.exists() {
        std::fs::remove_dir_all(&tmp)?;
    }

    // the decoys: every recruit the judge did not reject, listed in
    // decoys.tbl. each tool searches only its gaps, the decoys it did not
    // recruit itself, laid out per tool and family the way gather laid out
    // the recruits; learn takes the rest of its scores from recruit. the
    // list is streamed to disk rather than held: this process is about to
    // launch 40,000 commands, and every launch forks it, page tables and all,
    // so what it holds in memory is what every launch costs
    let decoys_dir = layout.decoys();
    if decoys_dir.exists() {
        std::fs::remove_dir_all(&decoys_dir)?;
    }
    for tool in GAPS {
        std::fs::create_dir_all(decoys_dir.join(tool))?;
    }
    let mut recruits = RecruitsByFamily::open(&layout.recruits_tbl())?;

    let decoys_tbl = layout.decoys_tbl();
    let schema = toil::Schema::new(["family", "target"]);
    let mut widths = schema.widths();
    widths.widen(&schema.row([
        "x".repeat(all_families.iter().map(String::len).max().unwrap_or(0)),
        "x".repeat(NAME_WIDTH),
    ]));
    let mut listed = toil::Stream::new(
        schema,
        widths,
        std::io::BufWriter::with_capacity(
            1 << 20,
            std::fs::File::create(&decoys_tbl)
                .with_context(|| format!("failed to create {}", decoys_tbl.display()))?,
        ),
    );

    // the families with a gap, per tool in GAPS order
    let mut gaps: [Vec<String>; 2] = Default::default();
    let (mut kept, mut dropped, mut families) = (0usize, 0usize, 0usize);
    for family in all_families {
        let tools = recruits.family(&family)?;
        let mut outs: [Option<std::io::BufWriter<std::fs::File>>; 2] = Default::default();
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
            let bits = *tools.get(&name).with_context(|| {
                format!("{name} is a recruit of {family} that recruits.tbl does not list")
            })?;
            listed.row([family.as_str(), name.as_str()])?;
            n += 1;

            for (out, tool) in outs.iter_mut().zip(GAPS) {
                if bits & bit(tool) != 0 {
                    continue;
                }
                let out = match out {
                    Some(out) => out,
                    None => out.insert(std::io::BufWriter::new(std::fs::File::create(
                        decoys_dir.join(tool).join(format!("{family}.fa")),
                    )?)),
                };
                out.write_all(rec)?;
                if !rec.ends_with(b"\n") {
                    out.write_all(b"\n")?;
                }
            }
        }
        kept += n;
        families += usize::from(n > 0);

        let mut gapped = false;
        for (i, out) in outs.into_iter().enumerate() {
            if let Some(mut out) = out {
                out.flush()?;
                gaps[i].push(family.clone());
                gapped = true;
            }
        }
        if gapped {
            // made here, in this thread, rather than as a command: a mkdir
            // is not worth a fork of this process
            std::fs::create_dir_all(tmp.join(&family))?;
        }
    }
    listed.header()?;
    drop(listed);
    drop(rejects);
    drop(recruits);
    println!("{kept} decoys over {families} families, {dropped} recruits rejected");
    for (tool, families) in GAPS.iter().zip(&gaps) {
        println!("{tool} fills gaps in {} families", families.len());
    }
    if families == 0 {
        bail!("every recruit was rejected; nothing to fill");
    }

    let nail_bin = nail()?;
    let mmseqs_bin = mmseqs()?;
    let jobs = args.jobs;
    let queries = layout.queries();
    let target = |tool: &str, family: &str| decoys_dir.join(tool).join(format!("{family}.fa"));
    let hmm = |family: &str| queries.join(family).join("query.hmm");
    let table = |tool: &str, family: &str| manifest::table_path(&results, tool, family);
    let scratch = |family: &str| tmp.join(family);

    println!("filling {jobs} families at once...");

    // a launch per family per tool with a gap there. mmseqs' chain -- its profile
    // from the family's alignment, its target db, the search, the table --
    // is five of its own processes, run as one shell command so that michi
    // launches it once and the manifest times it as one run
    let mmseqs_chain = |family: &str| -> String {
        let d = scratch(family);
        let m = mmseqs_bin.display();
        let q = |rest: &str| format!("{}/{rest}", d.display());
        [
            format!(
                "{m} convertmsa {} {} --identifier-field 0",
                queries.join(family).join("query.sto").display(),
                q("msaDB")
            ),
            format!(
                "{m} msa2profile --threads 1 {} {} --match-mode 1",
                q("msaDB"),
                q("queryDB")
            ),
            format!(
                "{m} createdb --threads 1 {} {}",
                target(MMSEQS, family).display(),
                q("targetDB")
            ),
            format!(
                "{m} search --threads 1 -k {} -s {FILL_S} --max-seqs {FILL_MAX_SEQS} -e {FILL_E} {} {} {} {}",
                util::search::MMSEQS_K,
                q("queryDB"),
                q("targetDB"),
                q("alnDB"),
                q("work")
            ),
            format!(
                "{m} convertalis --threads 1 --format-mode 0 {} {} {} {}",
                q("queryDB"),
                q("targetDB"),
                q("alnDB"),
                table(MMSEQS, family).display()
            ),
        ]
        .join(" > /dev/null && ")
    };

    let [nail_gaps, mmseqs_gaps] = &gaps;
    let mut pl = PipelineBuilder::new().pool(jobs);
    if !nail_gaps.is_empty() {
        pl = pl.step(per_family(NAIL, jobs, nail_gaps, |family| {
            PCmd::new(&nail_bin)
                .sub("search")
                .arg("--mmseqs-path", &mmseqs_bin)
                .arg("-t", 1)
                .arg("--tmp-dir", scratch(family).join("nail"))
                .arg("--mmseqs-s", FILL_S)
                .arg("--seed-mode", FILL_SEED_MODE)
                .arg("-S", FILL_SEED_P)
                .arg("--mmseqs-max-seqs", FILL_MAX_SEQS)
                .arg("-E", FILL_E)
                .flag("--allow-overwrite")
                .arg("--tbl-out", table(NAIL, family))
                .path(hmm(family))
                .path(target(NAIL, family))
                .field(manifest::NAME, NAIL)
                .field(manifest::TOOL, NAIL)
                .field(manifest::SHARD, family)
        }));
    }
    if !mmseqs_gaps.is_empty() {
        pl = pl.step(per_family(MMSEQS, jobs, mmseqs_gaps, |family| {
            PCmd::new("/bin/sh")
                .arg("-c", mmseqs_chain(family))
                .field(manifest::NAME, MMSEQS)
                .field(manifest::TOOL, MMSEQS)
                .field(manifest::SHARD, family)
        }));
    }
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
    ledger::record(&stage.root)?;

    // one rm for every family's scratch, in this thread, for the same reason
    // the mkdirs were. stderr stays: it is where a failed command's output is
    for family in nail_gaps.iter().chain(mmseqs_gaps) {
        std::fs::remove_dir_all(scratch(family)).ok();
    }
    Ok(())
}

/// The tools `fill` searches, each over the decoys it did not recruit.
//
// not hmmer: its recruit search scored every pair it could report over a
// whole shard, so a decoy it did not recruit is one it would not score
const GAPS: [&str; 2] = [NAIL, MMSEQS];

/// `recruits.tbl` read a family at a time, in the order gather wrote it.
struct RecruitsByFamily {
    rows: toil::Reader<std::io::BufReader<std::fs::File>>,
    /// Where each tool's column is.
    columns: [usize; 3],
    /// The row read and not yet handed out: family, target, tool bits.
    next: Option<(String, String, u8)>,
}

impl RecruitsByFamily {
    fn open(path: &Path) -> anyhow::Result<Self> {
        let rows = toil::Reader::new(std::io::BufReader::with_capacity(
            1 << 20,
            std::fs::File::open(path).with_context(|| {
                format!(
                    "no recruits at {}; run `cutoffs gather` first",
                    path.display()
                )
            })?,
        ))?;
        let columns = TOOLS.map(|tool| rows.header().index(tool));
        let columns = match columns {
            [Some(a), Some(b), Some(c)] => [a, b, c],
            _ => bail!("{} has no column per tool", path.display()),
        };
        let mut recruits = RecruitsByFamily {
            rows,
            columns,
            next: None,
        };
        recruits.advance()?;
        Ok(recruits)
    }

    fn advance(&mut self) -> anyhow::Result<()> {
        self.next = None;
        while let Some((_, entry)) = self.rows.next_entry()? {
            let toil::Entry::Row(cells) = entry else {
                continue;
            };
            let text =
                |i: usize| String::from_utf8_lossy(cells.field(i).unwrap_or(b"")).into_owned();
            let mut bits = 0;
            for (tool, at) in TOOLS.iter().zip(self.columns) {
                if cells.field(at) == Some(b"1") {
                    bits |= bit(tool);
                }
            }
            self.next = Some((text(0), text(1), bits));
            return Ok(());
        }
        Ok(())
    }

    /// The recruits of `family` and the tools that reported each, asked for
    /// in sorted order.
    fn family(&mut self, family: &str) -> anyhow::Result<HashMap<String, u8>> {
        let mut out = HashMap::new();
        while let Some((f, target, bits)) = &self.next {
            match f.as_str().cmp(family) {
                std::cmp::Ordering::Less => {}
                std::cmp::Ordering::Equal => {
                    out.insert(target.clone(), *bits);
                }
                std::cmp::Ordering::Greater => break,
            }
            self.advance()?;
        }
        Ok(out)
    }
}

// ------------------------------------------------------------------- learn

/// How many decoy scores are recorded per family per tool.
const N_SCORES: usize = 5;

/// The tools a calibration scores, in the order `cutoffs.tbl` writes them.
const TOOLS: [&str; 3] = [NAIL, MMSEQS, HMMER];

/// The column block of hmmer's cutoffs learned from its best domain score.
const HMMER_DOM: &str = "hmmer_dom";

/// hmmer's `--tblout` read at its best single domain's score and E-value.
//
// the same rows as HmmerTable, so one fill table learns both
// sets of cutoffs: the full sequence score at column 5 and
// the best domain's at column 8
struct HmmerBestDomain;

impl HitColumns for HmmerBestDomain {
    const QUERY: usize = HmmerTable::QUERY;
    const TARGET: usize = HmmerTable::TARGET;
    const Q_START: Option<usize> = None;
    const Q_END: Option<usize> = None;
    const T_START: Option<usize> = None;
    const T_END: Option<usize> = None;
    const SCORE: usize = 8;
    const E_VALUE: usize = 7;
    const N_COLUMNS: usize = HmmerTable::N_COLUMNS;
}

/// One column of a family's null: its best decoy scores, highest first, and
/// how many decoys it scored.
#[derive(Default)]
struct Null {
    top: Vec<f32>,
    n: usize,
}

impl Null {
    fn add(&mut self, score: f32) {
        self.n += 1;
        self.keep(score);
    }

    fn keep(&mut self, score: f32) {
        let at = self.top.partition_point(|s| *s >= score);
        if at < N_SCORES {
            self.top.insert(at, score);
            self.top.truncate(N_SCORES);
        }
    }

    fn merge(&mut self, other: Null) {
        self.n += other.n;
        for s in other.top {
            self.keep(s);
        }
    }

    /// The scores padded to `N_SCORES`, and the count.
    //
    // a family with fewer than N_SCORES decoys pads with zero, which the
    // readers of cutoffs.tbl take as "no usable cutoff"
    fn cells(&self) -> (Vec<f32>, usize) {
        let scores = self
            .top
            .iter()
            .copied()
            .chain(std::iter::repeat(0.0))
            .take(N_SCORES)
            .collect();
        (scores, self.n)
    }
}

/// A family's null in each column of `cutoffs.tbl`: the tools in `TOOLS`
/// order, then hmmer's best domain.
type Nulls = [Null; 4];

/// Where hmmer's best-domain null sits in `Nulls`.
const DOM: usize = 3;

/// Each pair's best score in one hit table.
//
// best rather than last: a tool can report a pair more than once, and a
// decoy's score is the highest it reached
fn best_per_pair<T: HitColumns>(path: &Path) -> anyhow::Result<HashMap<(String, String), f32>> {
    let table = Table::<HitParser<T>>::open(path)
        .with_context(|| format!("failed to read {}", path.display()))?;
    let mut best: HashMap<(String, String), f32> = HashMap::new();
    for hit in table.iter() {
        let at = best
            .entry((hit.query.clone(), hit.target.clone()))
            .or_insert(hit.score);
        *at = at.max(hit.score);
    }
    Ok(best)
}

/// Fold one table's pairs into the nulls, leaving out the rejects.
fn fold<T: HitColumns>(
    path: &Path,
    column: usize,
    rejects: &HashSet<u128>,
    nulls: &mut HashMap<String, Nulls>,
) -> anyhow::Result<()> {
    for ((family, target), score) in best_per_pair::<T>(path)? {
        if rejects.contains(&pair_key(family.as_bytes(), target.as_bytes())) {
            continue;
        }
        nulls.entry(family).or_default()[column].add(score);
    }
    Ok(())
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
    let filled: Vec<(String, String)> = searched
        .runs()
        .filter_map(|row| {
            Some((
                row.get(manifest::NAME)?.to_string(),
                row.get(manifest::SHARD)?.to_string(),
            ))
        })
        .collect();

    let rejects: HashSet<u128> = read_rejects(&layout.rejects_tbl())?
        .iter()
        .map(|(f, t)| pair_key(f.as_bytes(), t.as_bytes()))
        .collect();

    // a decoy's score comes from recruit wherever the tool recruited it and
    // from fill's search of its gaps otherwise; each pair is in one shard,
    // and a gap is a pair the tool's recruit table does not hold, so no pair
    // is counted twice. ranked by score, since the two stages' E-values are
    // on different scales
    let recruit = layout.recruit()?.results();
    let pool = pool(args.threads)?;
    let mut nulls: HashMap<String, Nulls> = pool.install(|| {
        layout
            .targets()
            .par_iter()
            .map(|(shard, _)| -> anyhow::Result<HashMap<String, Nulls>> {
                let mut nulls = HashMap::new();
                let path = |tool: &str| manifest::table_path(&recruit, tool, shard);
                fold::<NailTable>(&path(NAIL), 0, &rejects, &mut nulls)?;
                fold::<BlastTable>(&path(MMSEQS), 1, &rejects, &mut nulls)?;
                fold::<HmmerTable>(&path(HMMER), 2, &rejects, &mut nulls)?;
                fold::<HmmerBestDomain>(&path(HMMER), DOM, &rejects, &mut nulls)?;
                Ok(nulls)
            })
            .try_reduce(HashMap::new, |mut a, b| {
                for (family, nulls) in b {
                    let into = a.entry(family).or_default();
                    for (into, null) in into.iter_mut().zip(nulls) {
                        into.merge(null);
                    }
                }
                Ok(a)
            })
    })?;

    for (tool, family) in &filled {
        let path = manifest::table_path(&results, tool, family);
        match tool.as_str() {
            NAIL => fold::<NailTable>(&path, 0, &rejects, &mut nulls)?,
            MMSEQS => fold::<BlastTable>(&path, 1, &rejects, &mut nulls)?,
            _ => bail!("fill ran {tool:?}, which only fills nail's and mmseqs' gaps"),
        }
    }

    let mut headers = vec!["family".to_string()];
    for tool in TOOLS {
        headers.extend((1..=N_SCORES).map(|i| format!("{tool}_{i}")));
        headers.push(format!("{tool}_n"));
    }
    headers.extend((1..=N_SCORES).map(|i| format!("{HMMER_DOM}_{i}")));

    let mut families: Vec<(String, Nulls)> = nulls.into_iter().collect();
    families.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    let mut table = toil::Table::new(toil::Schema::new(headers));
    for (family, nulls) in families {
        let mut row = vec![family];
        for null in &nulls[..DOM] {
            row.extend(cells(&null.cells()));
        }
        // the same decoys scored by their best domain, so the count is the
        // one hmmer_n already holds
        row.extend(cells(&nulls[DOM].cells()).into_iter().take(N_SCORES));
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
    let jobs = args.jobs.context("--jobs is required")?;

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
            threads,
        },
        paths,
    )?;

    reject(
        RejectArgs {
            place: args.place.clone(),
            dry_run: false,
            jobs,
            reject_e: REJECT_E,
            reject_with: None,
            reject_set: None,
        },
        paths,
    )?;

    fill(
        FillArgs {
            place: args.place.clone(),
            dry_run: false,
            jobs,
        },
        paths,
    )?;

    learn(
        LearnArgs {
            place: args.place.clone(),
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
