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

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use anyhow::{bail, Context};
use clap::{Parser, Subcommand};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use bioio::tbl::{BlastTable, HitTable, HmmerTable, NailTable};
use run::{Asset, Bin, Ctx, Numa, Options, Search};

use crate::build;

/// Name of the calibration directory when none is given.
pub const DEFAULT_OUT: &str = "cutoffs";

/// Config blocks for the two search stages, both in `cutoffs.toml`.
const RECRUIT_BLOCK: &str = "recruit";
const DECOY_BLOCK: &str = "decoy";

// ------------------------------------------------------------------ layout

/// Where every artifact of a calibration lives.
///
/// The whole thing hangs off one directory inside the benchmark, so several
/// calibrations — different shard counts, different parameters — can sit side
/// by side and be deleted as a unit.
struct Layout {
    bench: PathBuf,
    root: PathBuf,
}

impl Layout {
    fn new(name: &str, out: &str) -> anyhow::Result<Self> {
        let bench = build::dir().join(name);
        if !bench.is_dir() {
            bail!(
                "benchmark directory {} does not exist; run `mgnify build{}` first",
                bench.display(),
                if name == build::DEFAULT_NAME {
                    String::new()
                } else {
                    format!(" --name {name}")
                }
            );
        }

        let bench = bench.canonicalize()?;
        let root = bench.join(out);
        Ok(Layout { bench, root })
    }

    /// Forward shards, built by `mgnify build`.
    fn mgy(&self) -> PathBuf {
        self.bench.join("mgy")
    }

    fn query_hmm(&self) -> PathBuf {
        self.bench.join("query.hmm")
    }

    fn query_sto(&self) -> PathBuf {
        self.bench.join("query.sto")
    }

    fn mgy_rev(&self) -> PathBuf {
        self.root.join("mgy-rev")
    }

    fn recruit(&self) -> PathBuf {
        self.root.join("recruit")
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
    /// Reverse the benchmark's shards into the calibration directory.
    Reverse(ReverseArgs),
    /// Sweep every family against the reversed shards to find candidates.
    Recruit(StageArgs),
    /// Un-reverse what hit, group it per family, and split the query set.
    Decoys(DecoysArgs),
    /// Search each family against its own decoys, forward and reversed.
    Search(StageArgs),
    /// Turn the decoy scores into per-family cutoffs.
    Learn(LearnArgs),
    /// Run every stage in order.
    All(AllArgs),
}

/// Which benchmark, and which calibration directory inside it.
#[derive(Parser, Debug, Clone)]
pub struct Where {
    /// Which benchmark directory under benchmarks/mgnify/ to read.
    #[arg(long, default_value = build::DEFAULT_NAME)]
    pub name: String,

    /// Calibration directory, created under the benchmark.
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

/// Arguments shared by the two stages that actually run search tools.
#[derive(Parser, Debug)]
pub struct StageArgs {
    #[command(flatten)]
    pub place: Where,

    /// Only run entries whose name matches this glob.
    #[arg(short, long)]
    pub filter: Option<String>,

    /// Override the thread count from cutoffs.toml.
    #[arg(short, long)]
    pub threads: Option<usize>,

    /// How many searches to keep in flight at once.
    #[arg(short = 'j', long)]
    pub jobs: Option<usize>,

    /// Pin to a NUMA node. Absent means no pinning and no numactl call.
    #[arg(long)]
    pub numa_node: Option<usize>,

    /// List the expanded runs and exit without executing anything.
    #[arg(long)]
    pub dry_run: bool,
}

#[derive(Parser, Debug)]
pub struct DecoysArgs {
    #[command(flatten)]
    pub place: Where,

    /// Run name to read for each tool. Only needed when the recruit config
    /// swept a tool over several settings.
    #[arg(long, value_name = "RUN")]
    pub nail: Option<String>,

    #[arg(long, value_name = "RUN")]
    pub mmseqs: Option<String>,

    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Parser, Debug)]
pub struct LearnArgs {
    #[command(flatten)]
    pub place: Where,

    #[arg(long, value_name = "RUN")]
    pub nail: Option<String>,

    #[arg(long, value_name = "RUN")]
    pub mmseqs: Option<String>,

    #[arg(long, value_name = "RUN")]
    pub hmmer: Option<String>,

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

    #[arg(long)]
    pub numa_node: Option<usize>,
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

/// Shard files named `<n>.fa`, in numeric order.
///
/// The index has to come off the stem as a whole: run names embed floating
/// point parameters, so splitting on dots is not safe elsewhere and is not
/// worth doing differently here.
fn shards(dir: &Path) -> anyhow::Result<Vec<(usize, PathBuf)>> {
    let entries = std::fs::read_dir(dir)
        .with_context(|| format!("failed to read {}", dir.display()))?;

    let mut out: Vec<(usize, PathBuf)> = entries
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "fa"))
        .filter_map(|p| {
            let i = p.file_stem()?.to_str()?.parse::<usize>().ok()?;
            Some((i, p))
        })
        .collect();

    out.sort_by_key(|(i, _)| *i);

    if out.is_empty() {
        bail!("no shards named <n>.fa in {}", dir.display());
    }

    Ok(out)
}

fn reverse(args: ReverseArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.name, &args.place.out)?;
    let src = layout.mgy();
    let dst = layout.mgy_rev();

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
        found.par_iter().try_for_each(|(i, path)| -> anyhow::Result<()> {
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

fn recruit(args: StageArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.name, &args.place.out)?;
    let rev = layout.mgy_rev();

    if !rev.is_dir() {
        bail!(
            "no reversed shards in {}; run `mgnify cutoffs reverse` first",
            rev.display()
        );
    }

    let searches: Vec<Search> = shards(&rev)?
        .into_iter()
        .map(|(i, path)| {
            Search::new(i.to_string(), path)
                .with(Asset::Hmm, layout.query_hmm())
                .with(Asset::Sto, layout.query_sto())
        })
        .collect();

    // one search per shard is already large enough to use every core
    let into = layout.recruit();
    stage(RECRUIT_BLOCK, searches, &args, into, 1)
}

// ------------------------------------------------------------------ decoys

fn decoys(args: DecoysArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.name, &args.place.out)?;
    let recruit_dir = layout.recruit();

    let runs = run::Runs::load(
        recruit_dir.join(run::table::FILE_NAME),
        recruit_dir.join("results"),
    )
    .context("could not read the recruit stage's runs table")?;

    let nail = match &args.nail {
        Some(n) => n.as_str(),
        None => runs.only_for_tool("nail")?,
    };
    let mmseqs = match &args.mmseqs {
        Some(n) => n.as_str(),
        None => runs.only_for_tool("mmseqs")?,
    };

    let targets = runs.shared_targets(&[nail, mmseqs]);
    if targets.is_empty() {
        bail!("nail and mmseqs have no successfully recruited shards in common");
    }

    println!("reading {} recruited shards...", targets.len());

    let pool = pool(args.threads)?;

    // per shard, which families each recruited sequence hit. Only names travel
    // here; the sequences themselves are read in the second pass.
    let wanted: Vec<(String, HashMap<String, Vec<String>>)> = pool.install(|| {
        targets
            .par_iter()
            .map(|target| -> anyhow::Result<(String, HashMap<String, Vec<String>>)> {
                let mut map: HashMap<String, Vec<String>> = HashMap::new();

                let path = runs.table_path(nail, target);
                let tbl = HitTable::from_path::<_, NailTable>(&path)
                    .with_context(|| format!("failed to read {}", path.display()))?;
                collect(tbl, &mut map);

                let path = runs.table_path(mmseqs, target);
                let tbl = HitTable::from_path::<_, BlastTable>(&path)
                    .with_context(|| format!("failed to read {}", path.display()))?;
                collect(tbl, &mut map);

                for v in map.values_mut() {
                    v.sort();
                    v.dedup();
                }

                Ok((target.clone(), map))
            })
            .collect::<anyhow::Result<Vec<_>>>()
    })?;

    let families: HashSet<String> = wanted
        .iter()
        .flat_map(|(_, m)| m.keys().cloned())
        .collect();

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

    let mgy = layout.mgy();
    pool.install(|| {
        wanted.par_iter().try_for_each(|(target, map)| -> anyhow::Result<()> {
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

            let path = mgy.join(target);
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

fn search(args: StageArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.name, &args.place.out)?;
    let decoy_dir = layout.decoys();

    if !decoy_dir.is_dir() {
        bail!(
            "no decoys in {}; run `mgnify cutoffs decoys` first",
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

    // two searches per family. The `.rev` suffix is deliberate: outputs are
    // keyed by the target's stem, so these land as `<run>.<family>.tbl` and
    // `<run>.<family>.rev.tbl`, which is the pairing `learn` reads back.
    let mut searches = Vec::with_capacity(families.len() * 2);
    for family in &families {
        let hmm = hmm_dir.join(format!("{family}.hmm"));
        let sto = sto_dir.join(format!("{family}.sto"));

        searches.push(
            Search::new(family.clone(), decoy_dir.join(format!("{family}.fa")))
                .with(Asset::Hmm, &hmm)
                .with(Asset::Sto, &sto),
        );
        searches.push(
            Search::new(
                format!("{family}.rev"),
                rev_dir.join(format!("{family}.rev.fa")),
            )
            .with(Asset::Hmm, &hmm)
            .with(Asset::Sto, &sto),
        );
    }

    // each family is small, so the parallelism belongs on the outside
    let jobs = args
        .jobs
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(4, |n| n.get()));

    stage(DECOY_BLOCK, searches, &args, layout.root.clone(), jobs)
}

// ------------------------------------------------------- shared stage setup

/// Plan and execute one of the two search stages.
fn stage(
    block: &str,
    searches: Vec<Search>,
    args: &StageArgs,
    into: PathBuf,
    default_jobs: usize,
) -> anyhow::Result<()> {
    let config = run::Config::from_path_as(build::dir().join("cutoffs.toml"), block)?;

    let opts = Options {
        filter: args.filter.clone(),
        threads: args.threads,
        numa_node: args.numa_node,
        jobs: args.jobs.unwrap_or(default_jobs),
        dry_run: args.dry_run,
    };
    let runs = run::plan(&config, &opts)?;

    if opts.dry_run {
        run::describe(&runs, &searches);
        return Ok(());
    }

    let results = into.join("results");
    if results.exists() {
        std::fs::remove_dir_all(&results)?;
    }
    std::fs::create_dir_all(&results)?;

    let threads = runs.iter().map(|r| r.threads).max().unwrap_or(1);
    let numa = match args.numa_node.or(config.defaults.numa_node) {
        Some(node) => Some(Numa::new(node, threads * opts.jobs)?),
        None => None,
    };

    let ctx = Ctx {
        bin: Bin::new(build::repo().join("tools/bin")),
        tmp: into.join("tmp"),
        results,
        // above results/, so it survives the wipe at the start of a stage
        runs_table: into.join(run::table::FILE_NAME),
        numa,
    };

    run::execute(&config, &runs, &searches, &ctx, opts.jobs)
}

// ------------------------------------------------------------------- learn

/// How many decoy scores are recorded per family per tool.
const N_SCORES: usize = 5;

fn learn(args: LearnArgs) -> anyhow::Result<()> {
    let layout = Layout::new(&args.place.name, &args.place.out)?;
    let results = layout.results();

    let runs = run::Runs::load(layout.root.join(run::table::FILE_NAME), &results)
        .context("could not read the search stage's runs table")?;

    let nail = match &args.nail {
        Some(n) => n.clone(),
        None => runs.only_for_tool("nail")?.to_string(),
    };
    let mmseqs = match &args.mmseqs {
        Some(n) => n.clone(),
        None => runs.only_for_tool("mmseqs")?.to_string(),
    };
    // hmmer is optional: it is here to be compared against nail's numbers, not
    // because anything downstream needs it
    let hmmer = match &args.hmmer {
        Some(n) => Some(n.clone()),
        None => runs.only_for_tool("hmmer").ok().map(str::to_string),
    };

    let mut families: Vec<String> = std::fs::read_dir(layout.decoys())?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "fa"))
        .filter_map(|p| Some(p.file_stem()?.to_str()?.to_string()))
        .collect();
    families.sort();

    let out_path = layout.cutoffs_txt();
    let out = Mutex::new(BufWriter::new(File::create(&out_path)?));
    let skipped = std::sync::atomic::AtomicUsize::new(0);

    let pool = pool(args.threads)?;
    pool.install(|| {
        families.par_iter().try_for_each(|family| -> anyhow::Result<()> {
            let nail_scores = decoy_scores::<NailTable>(&results, &nail, family, args.reverse_e_cutoff)?;
            let mmseqs_scores =
                decoy_scores::<BlastTable>(&results, &mmseqs, family, args.reverse_e_cutoff)?;
            let hmmer_scores = match &hmmer {
                Some(run) => decoy_scores::<HmmerTable>(&results, run, family, args.reverse_e_cutoff)?,
                None => None,
            };

            let (Some(n), Some(m)) = (nail_scores, mmseqs_scores) else {
                // a family without both tables tells us nothing comparative
                skipped.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                return Ok(());
            };

            let mut line = format!("{family},{},{}", group("nail", &n), group("mmseqs", &m));
            if let Some(h) = hmmer_scores {
                line.push(',');
                line.push_str(&group("hmmer", &h));
            }

            writeln!(
                out.lock().expect("output mutex poisoned"),
                "{line}"
            )?;
            Ok(())
        })
    })?;

    out.into_inner()
        .expect("output mutex poisoned")
        .flush()?;

    let skipped = skipped.load(std::sync::atomic::Ordering::Relaxed);
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
    let stage_args = |jobs: Option<usize>| StageArgs {
        place: args.place.clone(),
        filter: None,
        threads: None,
        jobs,
        numa_node: args.numa_node,
        dry_run: false,
    };

    reverse(ReverseArgs {
        place: args.place.clone(),
        shards: args.shards,
        threads: args.threads,
    })?;

    recruit(stage_args(Some(1)))?;

    decoys(DecoysArgs {
        place: args.place.clone(),
        nail: None,
        mmseqs: None,
        threads: args.threads,
    })?;

    search(stage_args(args.jobs))?;

    learn(LearnArgs {
        place: args.place.clone(),
        nail: None,
        mmseqs: None,
        hmmer: None,
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
