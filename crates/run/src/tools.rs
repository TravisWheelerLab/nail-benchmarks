use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::config::Run;
use crate::exec::{self, Job, Numa, Timing};
use crate::split::{self, Kind};

/// A query artifact a benchmark can offer. Tools ask for what they need and
/// fail cleanly when a benchmark does not provide it — mgnify has no unaligned
/// query fasta, long-seqs has no profiles.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum Asset {
    /// HMMER3 profiles.
    Hmm,
    /// Stockholm alignments.
    Sto,
    /// Unaligned query sequences.
    Fasta,
    /// Directory of per-family aligned fasta.
    Afa,
}

impl Asset {
    fn describe(&self) -> &'static str {
        match self {
            Asset::Hmm => "an hmm profile file",
            Asset::Sto => "a stockholm alignment file",
            Asset::Fasta => "a query fasta file",
            Asset::Afa => "a directory of aligned fasta",
        }
    }
}

/// One unit of work: a set of query artifacts searched against one target.
///
/// The benchmark builds this list, so the runner never has to guess at a
/// directory layout. pct-id yields one, mgnify one per shard, and long-seqs
/// six zipped query/target pairs.
#[derive(Clone, Debug)]
pub struct Search {
    /// Distinguishes outputs when a benchmark has more than one search;
    /// empty when there is only one.
    pub label: String,
    pub target: PathBuf,
    assets: BTreeMap<Asset, PathBuf>,
}

impl Search {
    pub fn new(label: impl Into<String>, target: impl Into<PathBuf>) -> Self {
        Search {
            label: label.into(),
            target: target.into(),
            assets: BTreeMap::new(),
        }
    }

    pub fn with(mut self, asset: Asset, path: impl Into<PathBuf>) -> Self {
        self.assets.insert(asset, path.into());
        self
    }

    pub fn asset(&self, asset: Asset) -> Result<&Path> {
        self.assets
            .get(&asset)
            .map(PathBuf::as_path)
            .with_context(|| {
                format!(
                    "this benchmark provides no {}, which the requested tool and query mode need",
                    asset.describe()
                )
            })
    }

    /// How this search is identified in the runs table.
    pub fn display(&self) -> String {
        if self.label.is_empty() {
            self.target
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_else(|| self.target.display().to_string())
        } else {
            self.label.clone()
        }
    }
}

/// Resolved paths into tools/bin.
pub struct Bin {
    root: PathBuf,
}

impl Bin {
    /// Canonicalized so the `cmd` column of the runs table shows a path you can
    /// paste into a shell.
    pub fn new(root: impl Into<PathBuf>) -> Self {
        let root = root.into();
        let root = root.canonicalize().unwrap_or(root);
        Bin { root }
    }

    pub fn get(&self, name: &str) -> Result<PathBuf> {
        let path = self.root.join(name);
        if !path.exists() {
            bail!(
                "missing tool binary {}; run `make {}` from the repo root",
                path.display(),
                suggest_make_target(name)
            );
        }
        Ok(path)
    }
}

fn suggest_make_target(bin: &str) -> &'static str {
    match bin {
        "nail" => "nail",
        "hmmsearch" | "phmmer" | "hmmbuild" | "hmmemit" | "esl-seqstat" | "create-profmark" => {
            "hmmer"
        }
        "mmseqs" => "mmseqs",
        "blastp" | "psiblast" | "makeblastdb" => "blast",
        "lastal" | "lastdb" => "last",
        "diamond" => "diamond",
        _ => "setup",
    }
}

pub struct Ctx {
    pub bin: Bin,
    pub tmp: PathBuf,
    pub results: PathBuf,
    pub numa: Option<Numa>,
}

impl Ctx {
    fn numa(&self) -> Option<&Numa> {
        self.numa.as_ref()
    }

    pub fn log_dir(&self) -> PathBuf {
        self.results.join(".logs")
    }

    /// Output path for one run of one search.
    pub fn out(&self, run: &Run, search: &Search, ext: &str) -> PathBuf {
        if search.label.is_empty() {
            self.results.join(format!("{}.{ext}", run.name))
        } else {
            self.results
                .join(format!("{}.{}.{ext}", run.name, search.label))
        }
    }

    pub fn log_path(&self, run: &Run, search: &Search) -> PathBuf {
        let stem = if search.label.is_empty() {
            run.name.clone()
        } else {
            format!("{}.{}", run.name, search.label)
        };
        self.log_dir().join(format!("{stem}.err"))
    }

    /// A job for one benchmark run, with stderr captured to that run's log.
    fn job(&self, program: PathBuf, run: &Run, search: &Search) -> Job {
        Job::new(program).stderr_to(self.log_path(run, search))
    }

    /// A job for a setup step, attributed to the tool rather than to any run.
    fn prep_job(&self, program: PathBuf, tool: &str) -> Job {
        Job::new(program).stderr_to(self.log_dir().join(format!("prep-{tool}.err")))
    }

    /// Scratch directory scoped to a name, cleared before use.
    fn scratch(&self, name: &str) -> Result<PathBuf> {
        let dir = self.tmp.join(name);
        if dir.exists() {
            std::fs::remove_dir_all(&dir)?;
        }
        std::fs::create_dir_all(&dir)?;
        Ok(dir)
    }

    fn scratch_keep(&self, name: &str) -> Result<PathBuf> {
        let dir = self.tmp.join(name);
        std::fs::create_dir_all(&dir)?;
        Ok(dir)
    }
}

/// What a tool did for one run: how long it took and the command line used.
pub struct Outcome {
    pub timing: Timing,
    pub cmd: String,
}

pub trait Tool {
    /// Build databases or indices this tool needs for a given search. Called
    /// once per (tool, search); implementations skip work that already exists
    /// so a shared query set is not rebuilt for every shard.
    fn prep(&self, ctx: &Ctx, search: &Search) -> Result<()>;

    fn run(&self, ctx: &Ctx, search: &Search, run: &Run) -> Result<Outcome>;
}

pub fn get(name: &str) -> Result<Box<dyn Tool>> {
    Ok(match name {
        "nail" => Box::new(Nail),
        "hmmer" => Box::new(Hmmer),
        "mmseqs" => Box::new(Mmseqs),
        "blast" => Box::new(Blast),
        "last" => Box::new(Last),
        "diamond" => Box::new(Diamond),
        other => {
            bail!("unknown tool {other:?}; known tools: nail, hmmer, mmseqs, blast, last, diamond")
        }
    })
}

/// Which side of a tool to invoke. Every tool here has a profile mode and a
/// sequence mode; each maps `prf`/`seq` onto its own inputs and binaries.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Query {
    Profile,
    Sequence,
}

fn query_kind(run: &Run) -> Result<Query> {
    let raw = run.var_str("query").with_context(|| {
        format!(
            "run {:?} has no `query` key; every run block needs one to pick its input",
            run.name
        )
    })?;

    match raw.as_str() {
        "prf" => Ok(Query::Profile),
        "seq" => Ok(Query::Sequence),
        other => bail!(
            "query must be `prf` or `seq`, got {other:?} in run {:?}",
            run.name
        ),
    }
}

fn sequence_only(run: &Run, tool: &str) -> Result<()> {
    if query_kind(run)? != Query::Sequence {
        bail!(
            "{tool} has no profile mode; use query = \"seq\" in run {:?}",
            run.name
        );
    }
    Ok(())
}

/// A stable scratch name for a path, so databases derived from a shared query
/// set are built once rather than once per shard.
fn slug(path: &Path) -> String {
    path.to_string_lossy()
        .chars()
        .map(|c| if c.is_alphanumeric() { c } else { '_' })
        .collect::<String>()
        .trim_matches('_')
        .to_string()
}

// ---------------------------------------------------------------- nail

struct Nail;

impl Tool for Nail {
    fn prep(&self, _ctx: &Ctx, _search: &Search) -> Result<()> {
        Ok(())
    }

    fn run(&self, ctx: &Ctx, search: &Search, run: &Run) -> Result<Outcome> {
        let query = match query_kind(run)? {
            Query::Profile => search.asset(Asset::Hmm)?,
            Query::Sequence => search.asset(Asset::Fasta)?,
        };

        let tmp = ctx.scratch(&format!("nail/{}", run.name))?;
        let mmseqs = ctx.bin.get("mmseqs")?;

        let job = ctx
            .job(ctx.bin.get("nail")?, run, search)
            .arg("search")
            .arg("--mmseqs-path")
            .arg(mmseqs.display().to_string())
            .arg("-t")
            .arg(run.threads.to_string())
            .arg("--tmp-dir")
            .arg(tmp.display().to_string())
            .arg("--tbl-out")
            .arg(ctx.out(run, search, "tbl").display().to_string())
            .args(run.args.clone())
            .arg(query.display().to_string())
            .arg(search.target.display().to_string());

        let cmd = job.display(ctx.numa());
        let timing = exec::run(&job, ctx.numa())?;

        // nail drops seeds in its tmp dir; keep them beside the other outputs
        let seeds = tmp.join("seeds.tsv");
        if seeds.exists() {
            std::fs::rename(&seeds, ctx.out(run, search, "seeds"))?;
        }

        Ok(Outcome { timing, cmd })
    }
}

// --------------------------------------------------------------- hmmer

struct Hmmer;

impl Tool for Hmmer {
    fn prep(&self, _ctx: &Ctx, _search: &Search) -> Result<()> {
        Ok(())
    }

    fn run(&self, ctx: &Ctx, search: &Search, run: &Run) -> Result<Outcome> {
        let (program, query, split_kind) = match query_kind(run)? {
            Query::Profile => ("hmmsearch", search.asset(Asset::Hmm)?, Kind::Hmm),
            Query::Sequence => ("phmmer", search.asset(Asset::Fasta)?, Kind::Fasta),
        };

        let program = ctx.bin.get(program)?;
        let tbl = ctx.out(run, search, "tbl");

        // hmmsearch/phmmer scale poorly past a few threads, so the query set is
        // split and run as several concurrent processes when threads_per is set
        let splits = match run.threads_per {
            Some(per) if per < run.threads => run.threads / per,
            _ => 1,
        };

        if splits <= 1 {
            let job = ctx
                .job(program.clone(), run, search)
                .arg("--cpu")
                .arg(run.threads.to_string())
                .args(run.args.clone())
                .arg("-o")
                .arg("/dev/null")
                .arg("--tblout")
                .arg(tbl.display().to_string())
                .arg(query.display().to_string())
                .arg(search.target.display().to_string());

            let cmd = job.display(ctx.numa());
            let timing = exec::run(&job, ctx.numa())?;
            return Ok(Outcome { timing, cmd });
        }

        let tmp = ctx.scratch(&format!("hmmer/{}", run.name))?;
        let parts = split::write_splits(query, split_kind, splits, tmp.join("query"))?;
        let per = run.threads_per.expect("splits > 1 implies threads_per");

        let mut jobs = Vec::with_capacity(parts.len());
        let mut part_tbls = Vec::with_capacity(parts.len());

        for part in &parts {
            let part_tbl = part.with_extension("tbl");
            jobs.push(
                ctx.job(program.clone(), run, search)
                    .arg("--cpu")
                    .arg(per.to_string())
                    .args(run.args.clone())
                    .arg("-o")
                    .arg("/dev/null")
                    .arg("--tblout")
                    .arg(part_tbl.display().to_string())
                    .arg(part.display().to_string())
                    .arg(search.target.display().to_string()),
            );
            part_tbls.push(part_tbl);
        }

        let cmd = format!("[{} x] {}", jobs.len(), jobs[0].display(ctx.numa()));
        let timing = exec::run_concurrent(&jobs, ctx.numa())?;

        concat(&part_tbls, &tbl)?;

        Ok(Outcome { timing, cmd })
    }
}

// -------------------------------------------------------------- mmseqs

struct Mmseqs;

impl Mmseqs {
    fn target_db(ctx: &Ctx, search: &Search) -> PathBuf {
        ctx.tmp
            .join(format!("mmseqs/target-{}/db", slug(&search.target)))
    }

    fn query_db(ctx: &Ctx, search: &Search, query: Query) -> Result<PathBuf> {
        let (kind, source) = match query {
            Query::Profile => ("prf", search.asset(Asset::Sto)?),
            Query::Sequence => ("seq", search.asset(Asset::Fasta)?),
        };
        Ok(ctx
            .tmp
            .join(format!("mmseqs/query-{kind}-{}/db", slug(source))))
    }
}

impl Tool for Mmseqs {
    fn prep(&self, ctx: &Ctx, search: &Search) -> Result<()> {
        let mmseqs = ctx.bin.get("mmseqs")?;

        // target db is per-target, so it is rebuilt for every shard
        let tdb = Self::target_db(ctx, search);
        if !tdb.exists() {
            std::fs::create_dir_all(tdb.parent().expect("db path has a parent"))?;
            let job = ctx
                .prep_job(mmseqs.clone(), "mmseqs")
                .arg("createdb")
                .arg(search.target.display().to_string())
                .arg(tdb.display().to_string());
            check(&job, ctx.numa(), "mmseqs createdb (target)")?;
        }

        // query dbs are keyed by source path, so a query set shared across
        // shards is converted once rather than once per shard
        if let Ok(fasta) = search.asset(Asset::Fasta) {
            let qdb = Self::query_db(ctx, search, Query::Sequence)?;
            if !qdb.exists() {
                std::fs::create_dir_all(qdb.parent().expect("db path has a parent"))?;
                let job = ctx
                    .prep_job(mmseqs.clone(), "mmseqs")
                    .arg("createdb")
                    .arg(fasta.display().to_string())
                    .arg(qdb.display().to_string());
                check(&job, ctx.numa(), "mmseqs createdb (query)")?;
            }
        }

        if let Ok(sto) = search.asset(Asset::Sto) {
            let qdb = Self::query_db(ctx, search, Query::Profile)?;
            if !qdb.exists() {
                let dir = ctx.scratch_keep(&format!("mmseqs/query-prf-{}", slug(sto)))?;
                let msa_db = dir.join("msaDB");

                let job = ctx
                    .prep_job(mmseqs.clone(), "mmseqs")
                    .arg("convertmsa")
                    .arg(sto.display().to_string())
                    .arg(msa_db.display().to_string())
                    .arg("--identifier-field")
                    .arg("0");
                check(&job, ctx.numa(), "mmseqs convertmsa")?;

                let job = ctx
                    .prep_job(mmseqs.clone(), "mmseqs")
                    .arg("msa2profile")
                    .arg(msa_db.display().to_string())
                    .arg(qdb.display().to_string())
                    .arg("--match-mode")
                    .arg("1");
                check(&job, ctx.numa(), "mmseqs msa2profile")?;
            }
        }

        Ok(())
    }

    fn run(&self, ctx: &Ctx, search: &Search, run: &Run) -> Result<Outcome> {
        let mmseqs = ctx.bin.get("mmseqs")?;
        let qdb = Self::query_db(ctx, search, query_kind(run)?)?;
        let tdb = Self::target_db(ctx, search);

        // mmseqs aborts rather than overwrite an existing alignment db, so both
        // it and the working directory live in scratch that is cleared per
        // (run, search). Keying on the search too matters once a benchmark has
        // more than one, as mgnify's shards do.
        let key = if search.label.is_empty() {
            run.name.clone()
        } else {
            format!("{}.{}", run.name, search.label)
        };
        let adb = ctx.scratch(&format!("mmseqs/align/{key}"))?.join("db");
        let work = ctx.scratch(&format!("mmseqs/work/{key}"))?;

        let job = ctx
            .job(mmseqs.clone(), run, search)
            // mmseqs reports failures on stdout, not stderr, so it shares the
            // run log; the log is discarded when the run succeeds
            .stdout_append(ctx.log_path(run, search))
            .arg("search")
            .arg(qdb.display().to_string())
            .arg(tdb.display().to_string())
            .arg(adb.display().to_string())
            .arg(work.display().to_string())
            .arg("--threads")
            .arg(run.threads.to_string())
            .args(run.args.clone());

        let cmd = job.display(ctx.numa());
        let timing = exec::run(&job, ctx.numa())?;

        // convertalis is bookkeeping, not part of the measured search
        let job = ctx
            .job(mmseqs, run, search)
            .stdout_append(ctx.log_path(run, search))
            .arg("convertalis")
            .arg(qdb.display().to_string())
            .arg(tdb.display().to_string())
            .arg(adb.display().to_string())
            .arg(ctx.out(run, search, "tbl").display().to_string())
            .arg("--format-mode")
            .arg("0");
        check(&job, ctx.numa(), "mmseqs convertalis")?;

        std::fs::remove_dir_all(&work).ok();

        Ok(Outcome { timing, cmd })
    }
}

// --------------------------------------------------------------- blast

struct Blast;

impl Blast {
    fn db(ctx: &Ctx, search: &Search) -> PathBuf {
        ctx.tmp
            .join(format!("blast/{}/target_db", slug(&search.target)))
    }
}

impl Tool for Blast {
    fn prep(&self, ctx: &Ctx, search: &Search) -> Result<()> {
        let db = Self::db(ctx, search);
        if db.with_extension("pdb").exists() {
            return Ok(());
        }
        std::fs::create_dir_all(db.parent().expect("db path has a parent"))?;

        let job = ctx
            .prep_job(ctx.bin.get("makeblastdb")?, "blast")
            .arg("-in")
            .arg(search.target.display().to_string())
            .arg("-dbtype")
            .arg("prot")
            .arg("-out")
            .arg(db.display().to_string());
        check(&job, ctx.numa(), "makeblastdb")
    }

    fn run(&self, ctx: &Ctx, search: &Search, run: &Run) -> Result<Outcome> {
        let tbl = ctx.out(run, search, "tbl");

        match query_kind(run)? {
            Query::Sequence => {
                let job = ctx
                    .job(ctx.bin.get("blastp")?, run, search)
                    .arg("-query")
                    .arg(search.asset(Asset::Fasta)?.display().to_string())
                    .arg("-db")
                    .arg(Self::db(ctx, search).display().to_string())
                    .arg("-out")
                    .arg(tbl.display().to_string())
                    .arg("-outfmt")
                    .arg("6")
                    .arg("-num_threads")
                    .arg(run.threads.to_string())
                    .args(run.args.clone());

                let cmd = job.display(ctx.numa());
                let timing = exec::run(&job, ctx.numa())?;
                Ok(Outcome { timing, cmd })
            }
            // psiblast takes one MSA at a time, so a profile run is one
            // invocation per family with output appended to a single table
            Query::Profile => {
                let psiblast = ctx.bin.get("psiblast")?;
                let afa_dir = search.asset(Asset::Afa)?;

                let mut msas: Vec<PathBuf> = std::fs::read_dir(afa_dir)
                    .with_context(|| format!("failed to read {}", afa_dir.display()))?
                    .filter_map(|e| e.ok().map(|e| e.path()))
                    .filter(|p| p.extension().is_some_and(|e| e == "afa"))
                    .collect();
                msas.sort();

                if msas.is_empty() {
                    bail!("no .afa files in {}", afa_dir.display());
                }

                std::fs::remove_file(&tbl).ok();
                if let Some(dir) = tbl.parent() {
                    std::fs::create_dir_all(dir)?;
                }

                let mut parts = Vec::with_capacity(msas.len());
                let mut first_cmd = None;
                let start = std::time::Instant::now();

                for msa in &msas {
                    let job = ctx
                        .job(psiblast.clone(), run, search)
                        .arg("-in_msa")
                        .arg(msa.display().to_string())
                        .arg("-db")
                        .arg(Self::db(ctx, search).display().to_string())
                        .arg("-outfmt")
                        .arg("6")
                        .arg("-num_threads")
                        .arg(run.threads.to_string())
                        .arg("-comp_based_stats")
                        .arg("1")
                        .arg("-num_iterations")
                        .arg("1")
                        .args(run.args.clone())
                        .stdout_append(&tbl);

                    if first_cmd.is_none() {
                        first_cmd = Some(job.display(ctx.numa()));
                    }
                    parts.push(exec::run(&job, ctx.numa())?);
                }

                let timing = Timing::combine(&parts, start.elapsed().as_secs_f64());
                let cmd = format!("[{} x] {}", msas.len(), first_cmd.unwrap_or_default());
                Ok(Outcome { timing, cmd })
            }
        }
    }
}

// ---------------------------------------------------------------- last

struct Last;

impl Last {
    fn db(ctx: &Ctx, search: &Search) -> PathBuf {
        ctx.tmp
            .join(format!("last/{}/target_db", slug(&search.target)))
    }
}

impl Tool for Last {
    fn prep(&self, ctx: &Ctx, search: &Search) -> Result<()> {
        let db = Self::db(ctx, search);
        if db.with_extension("prj").exists() {
            return Ok(());
        }
        std::fs::create_dir_all(db.parent().expect("db path has a parent"))?;

        let job = ctx
            .prep_job(ctx.bin.get("lastdb")?, "last")
            .arg("-p")
            .arg(db.display().to_string())
            .arg(search.target.display().to_string());
        check(&job, ctx.numa(), "lastdb")
    }

    fn run(&self, ctx: &Ctx, search: &Search, run: &Run) -> Result<Outcome> {
        sequence_only(run, "last")?;

        // lastal writes its table to stdout
        let job = ctx
            .job(ctx.bin.get("lastal")?, run, search)
            .arg(Self::db(ctx, search).display().to_string())
            .arg(search.asset(Asset::Fasta)?.display().to_string())
            .arg("-f")
            .arg("BlastTab")
            .arg("-P")
            .arg(run.threads.to_string())
            .args(run.args.clone())
            .stdout_to(ctx.out(run, search, "tbl"));

        let cmd = job.display(ctx.numa());
        let timing = exec::run(&job, ctx.numa())?;
        Ok(Outcome { timing, cmd })
    }
}

// ------------------------------------------------------------- diamond

struct Diamond;

impl Diamond {
    fn db(ctx: &Ctx, search: &Search) -> PathBuf {
        ctx.tmp
            .join(format!("diamond/{}/target_db", slug(&search.target)))
    }
}

impl Tool for Diamond {
    fn prep(&self, ctx: &Ctx, search: &Search) -> Result<()> {
        let db = Self::db(ctx, search);
        if db.with_extension("dmnd").exists() {
            return Ok(());
        }
        std::fs::create_dir_all(db.parent().expect("db path has a parent"))?;

        let job = ctx
            .prep_job(ctx.bin.get("diamond")?, "diamond")
            .arg("makedb")
            .arg("--in")
            .arg(search.target.display().to_string())
            .arg("--db")
            .arg(db.display().to_string());
        check(&job, ctx.numa(), "diamond makedb")
    }

    fn run(&self, ctx: &Ctx, search: &Search, run: &Run) -> Result<Outcome> {
        sequence_only(run, "diamond")?;

        let job = ctx
            .job(ctx.bin.get("diamond")?, run, search)
            .arg("blastp")
            .arg("--query")
            .arg(search.asset(Asset::Fasta)?.display().to_string())
            .arg("--db")
            .arg(Self::db(ctx, search).display().to_string())
            .arg("--out")
            .arg(ctx.out(run, search, "tbl").display().to_string())
            .arg("--outfmt")
            .arg("6")
            .arg("--threads")
            .arg(run.threads.to_string())
            .args(run.args.clone());

        let cmd = job.display(ctx.numa());
        let timing = exec::run(&job, ctx.numa())?;
        Ok(Outcome { timing, cmd })
    }
}

// -------------------------------------------------------------- shared

/// Run a setup step and fail loudly if it does not succeed. Prep steps are not
/// part of the measured search, so their timing is discarded.
fn check(job: &Job, numa: Option<&Numa>, what: &str) -> Result<()> {
    let timing = exec::run(job, numa)?;
    if timing.exit != 0 {
        bail!(
            "{what} failed with exit code {}: {}",
            timing.exit,
            job.display(numa)
        );
    }
    Ok(())
}

fn concat(parts: &[PathBuf], out: &Path) -> Result<()> {
    if let Some(dir) = out.parent() {
        std::fs::create_dir_all(dir)?;
    }
    let mut dst =
        std::fs::File::create(out).with_context(|| format!("failed to create {}", out.display()))?;

    for part in parts {
        let mut src = std::fs::File::open(part)
            .with_context(|| format!("failed to open {}", part.display()))?;
        std::io::copy(&mut src, &mut dst)?;
    }

    Ok(())
}
