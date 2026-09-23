//! The steps every pipeline is assembled out of.
//!
//! These return [`Step`]s and [`Cmd`]s for a pipeline to compose. Nothing here
//! owns a pipeline or decides what a run measures: a pipeline searching hmmer
//! against a shard takes those steps and puts them where it needs them.
//!
//! Every search command carries the fields `parse` reads back: `name` for the
//! column it becomes, `tool` for how to read its table, and `shard` for which
//! target it ran against. A command that is not itself a search carries `stage`
//! instead -- seeding, which produces pairs rather than scores, and mmseqs'
//! conversion, which reformats a search's output. Only a `name` makes a
//! command a column, so only search time is charged to one.

use std::path::{Path, PathBuf};

use anyhow::Context;

use michi::{Closure, Cmd, Step};

use crate::manifest;
use crate::split::{self, Kind};

// the four settings below are held equal across the benchmarks so
// that their tools compare: a seed set from one pipeline is the
// seed set another would have got, and a hit one reports is one
// the others would have reported. what a benchmark sweeps is its
// own and lives with it

/// nail's seeding mode, for every pipeline that seeds.
pub const SEED_MODE: &str = "prog";

/// The sensitivity the pipelines that seed once seed at.
//
// cloud-search and loss-decomp used to spell this `MMSEQS_S` apart, with a
// comment on each saying it had to match the other
pub const SEED_S: &str = "12.0";

/// hmmsearch doesn't scale past a couple of threads, so its query gets split
/// threads/HMMER_CPU ways and the parts run at the same time.
pub const HMMER_CPU: usize = 2;

/// Every tool reports down to here, so they can be compared.
pub const EVALUE: &str = "10";

/// The stage the query split belongs to.
//
// no table and no search, but splitting all of Pfam takes long
// enough that a run omitting it does not account for its own
// wall clock
const SPLIT: &str = "split";

/// The stage building mmseqs' target database belongs to.
const CREATEDB: &str = "createdb";

/// The stage mmseqs' table conversion belongs to.
//
// a stage rather than a run: what a column cost is the search that produced
// the alignments, not the pass that reformats them. every tool is timed the
// same way -- one search command, nothing around it
const CONVERT: &str = "convert";

/// Where one pipeline's output lives.
///
/// Everything a run produces -- hit tables, domain tables, the seed list --
/// lands in `results/`, told apart by name rather than by directory. The
/// scratch it wanted on the way is a directory of its own, so what a run
/// produced and what it merely needed are never the same tree.
pub struct Dirs {
    pub root: PathBuf,
    results: PathBuf,
    pub tmp: PathBuf,
}

impl Dirs {
    /// Told where the run goes and where its scratch goes. This crate builds
    /// command lines; where their output lands is the caller's to know.
    pub fn new(run: impl Into<PathBuf>, tmp: impl Into<PathBuf>) -> Dirs {
        let root = run.into();
        Dirs {
            results: root.join("results"),
            tmp: tmp.into(),
            root,
        }
    }

    /// Makes everything a pipeline writes into. hmmsearch won't create its own
    /// output directory; it just fails to open its output.
    pub fn mkdir(&self) -> Cmd {
        Cmd::new("mkdir")
            .name("dirs")
            .flag("-p")
            .path(&self.results)
            .path(self.tmp.join("hmmer"))
    }

    /// Everything a run leaves behind, cleared at the front of the pipeline
    /// rather than before it.
    ///
    /// A stale table from an earlier run reads as a run that simply found
    /// less, and mmseqs refuses to overwrite an existing alignment db. Doing
    /// it as a step is what keeps `--dry-run` from touching the disk.
    pub fn clean(&self) -> Step {
        Step::serial([
            Cmd::new("rm")
                .name("clean")
                .flag("-rf")
                .path(&self.results)
                .path(&self.tmp),
            Cmd::new("mkdir")
                .name("dirs")
                .flag("-p")
                .path(&self.results),
        ])
        .name("clean")
    }

    pub fn table(&self, name: &str, shard: &str) -> PathBuf {
        manifest::table_path(&self.results, name, shard)
    }

    /// Where one seeding's list for one shard goes.
    pub fn seeds(&self, seeding: &str, shard: &str) -> PathBuf {
        manifest::seeds_path(&self.results, seeding, shard)
    }
}

/// What a step is called: the run, and the shard where there is one.
fn step_name(name: &str, shard: &str) -> String {
    match shard.is_empty() {
        true => name.to_string(),
        false => format!("{name}.{shard}"),
    }
}

/// Tag one command with what the ledger reads back, and wrap it as a step.
///
/// The fields every search carries: the run it belongs to and the tool that
/// produced it. `extra` is everything else the ledger reads back -- which
/// shard, which mode, what a sweep moved.
pub fn tag(cmd: Cmd, name: &str, tool: &str, extra: &[(&str, String)]) -> Cmd {
    extra.iter().fold(
        cmd.field(manifest::NAME, name).field(manifest::TOOL, tool),
        |cmd, (key, value)| cmd.field(*key, value),
    )
}

/// nail's search, as every benchmark here runs it.
///
/// `args` is what this benchmark is moving: a sensitivity, an `(A, B)` cell,
/// a seeding to replay. `flags` is the same for the ones that take no value.
/// `fields` is what the ledger reads back beside the run's name.
#[allow(clippy::too_many_arguments)]
pub fn nail(
    bin: &Path,
    mmseqs: &Path,
    query: &Path,
    target: &Path,
    out: &Path,
    tmp: &Path,
    threads: usize,
    evalue: &str,
    args: &[(&str, String)],
    flags: &[&str],
    name: &str,
    fields: &[(&str, String)],
) -> Step {
    let cmd = Cmd::new(bin)
        .sub("search")
        // nail looks for mmseqs at startup even where it will never call it,
        // and nothing here is on PATH
        .arg("--mmseqs-path", mmseqs)
        .arg("-t", threads)
        .arg("--tmp-dir", tmp)
        .flag("--allow-overwrite")
        .arg("-E", evalue)
        .arg("--tbl-out", out);

    let cmd = args
        .iter()
        .fold(cmd, |cmd, (key, value)| cmd.arg(*key, value));
    let cmd = flags.iter().fold(cmd, |cmd, flag| cmd.flag(*flag));

    Step::serial([tag(cmd.path(query).path(target), name, "nail", fields)]).name(name)
}

/// There's no shell to expand a glob, so the parts get named one by one.
fn cat(parts: impl IntoIterator<Item = PathBuf>, into: PathBuf) -> Cmd {
    parts
        .into_iter()
        .fold(Cmd::new("cat"), |cmd, part| cmd.path(part))
        .stdout_to(into)
}

/// A query set cut into parts for hmmer to search in parallel.
///
/// The parts are named before they exist: `write_splits` files them by index,
/// so a batch over them can be written without waiting to see what the split
/// produced. Splitting is separate from searching because the parts don't
/// depend on the target, and a pipeline that searches many shards splits once.
pub struct Split {
    query: PathBuf,
    kind: Kind,
    dir: PathBuf,
    parts: Vec<PathBuf>,
}

impl Split {
    pub fn new(
        query: impl Into<PathBuf>,
        kind: Kind,
        dir: impl Into<PathBuf>,
        jobs: usize,
    ) -> Split {
        let dir = dir.into();
        let ext = kind.extension();

        Split {
            query: query.into(),
            kind,
            parts: (0..jobs).map(|i| dir.join(format!("{i}.{ext}"))).collect(),
            dir,
        }
    }

    /// Rust in place of a command, so it is a closure step. Whatever a previous
    /// run left in there would be searched as if it belonged.
    ///
    /// `name` is the step's, which is what a ledger row falls back to when a
    /// command carries no run name. `extra` is whatever tells one split from
    /// another, for a pipeline that cuts up more than one query set.
    pub fn step(&self, name: &str, extra: &[(&str, String)]) -> Step {
        let (query, kind, dir, jobs) = (
            self.query.clone(),
            self.kind,
            self.dir.clone(),
            self.parts.len(),
        );

        let closure = Closure::new("split", move || {
            std::fs::remove_dir_all(&dir).ok();
            let written = split::write_splits(&query, kind, jobs, &dir)?;

            // empty bins are skipped, so a query with fewer models than parts
            // comes back short and the batch would be pointed at files that
            // were never written
            anyhow::ensure!(
                written.len() == jobs,
                "split {} into {} parts, expected {jobs}",
                query.display(),
                written.len()
            );

            Ok(())
        })
        .field(manifest::STAGE, SPLIT);

        let closure = extra
            .iter()
            .fold(closure, |closure, (key, value)| closure.field(*key, value));

        Step::from_closures([closure]).name(name)
    }
}

/// The mmseqs database a search reads its targets out of.
pub fn createdb(mmseqs: &Path, target: &Path, db: &Path, shard: &str, threads: usize) -> Cmd {
    Cmd::new(mmseqs)
        .name("createdb")
        .sub("createdb")
        // mmseqs takes every core it can find unless told otherwise, so the
        // setup around a search is held to the same count as the search itself
        .arg("--threads", threads)
        .path(target)
        .path(db)
        // a stage, not a run: charging database construction to
        // a search would report the tool as slower than it is,
        // and dropping it would lose wall clock the pipeline paid
        .field(manifest::STAGE, CREATEDB)
        .field(manifest::SHARD, shard)
}

/// One hmmer run over one shard: the query's parts searched together, then
/// their output gathered up.
///
/// Named rather than a pair, so a caller that relabels the steps can say which
/// one it is relabelling.
pub struct Hmmer {
    pub search: Step,
    pub cat: Step,
}

/// One hmmer run over one shard, and the two tables it leaves behind.
///
/// Two steps: the query's parts searched together, then their output
/// concatenated into `results/<name>.<shard>.tbl` and `.domtbl`. It is an
/// ordinary run and carries the ordinary fields -- the domain table is the
/// only thing about it the other tools have no equivalent of.
///
/// `extra` is whatever else tells this run apart from the pipeline's others --
/// a repetition index, a swept setting. Without it a pipeline that runs hmmer
/// more than once against the same shard writes rows it cannot tell apart.
#[allow(clippy::too_many_arguments)]
pub fn hmmer(
    program: &Path,
    split: &Split,
    dirs: &Dirs,
    name: &str,
    tool: &str,
    shard_name: &str,
    target: &Path,
    evalue: &str,
    dom: bool,
    extra: &[(&str, String)],
) -> Hmmer {
    let parts = &split.parts;

    // beside the parts rather than in a directory of its own, so two runs
    // splitting two different queries do not write over each other
    let scratch = split.dir.clone();

    let fields = |cmd: Cmd| {
        let cmd = tag(cmd, name, tool, extra);
        match shard_name.is_empty() {
            true => cmd,
            false => cmd.field(manifest::SHARD, shard_name),
        }
    };

    Hmmer {
        search: Step::batched(
            parts.len(),
            parts.iter().enumerate().map(|(i, part)| {
                fields(
                    Cmd::new(program)
                        .name(i.to_string())
                        .arg("--cpu", HMMER_CPU)
                        .arg("--tblout", scratch.join(format!("{i}.tbl")))
                        .arg("--domtblout", scratch.join(format!("{i}.domtbl")))
                        .arg("-E", evalue)
                        .path(part)
                        .path(target),
                )
            }),
        )
        .name(step_name(name, shard_name))
        // per command, not per step, so this asks for HMMER_CPU x parts, which
        // is --threads again. a machine with a smaller pool than that won't
        // fail, it will just run fewer of the parts at once
        .cores(HMMER_CPU),
        cat: Step::serial(
            [fields(
                cat(
                    (0..parts.len()).map(|i| scratch.join(format!("{i}.tbl"))),
                    manifest::table_path(&dirs.results, name, shard_name),
                )
                .name("tbl"),
            )]
            .into_iter()
            .chain(dom.then(|| {
                fields(
                    cat(
                        (0..parts.len()).map(|i| scratch.join(format!("{i}.domtbl"))),
                        manifest::dom_path(&dirs.results, name, shard_name),
                    )
                    .name("domtbl"),
                )
            })),
        )
        .name(format!("cat.{}", step_name(name, shard_name))),
    }
}

/// The k-mer length every mmseqs search here is held to.
//
// nail passes this to the mmseqs prefilter it seeds with --
// `--mmseqs-k`, whose default is 6 -- while mmseqs' own
// default is 0, meaning it picks one from the database size.
// a standalone run left on 0 would be searching with a
// different k at every rung of the ladder, and a different
// one again from the k inside nail
const MMSEQS_K: usize = 6;

/// The stage a seeding belongs to, for the pipelines that record it as one.
//
// a ledger row is keyed by (stage, shard), so a pipeline that
// seeds several times per shard has to suffix this
pub const SEED: &str = "seed";

/// One seeding pass, kept so later searches can replay it.
///
/// The seeds are per-pair, so `parse` reads them back as the `seeded` column —
/// which is what lets a pipeline ask where a hit was lost rather than only
/// whether it survived.
#[allow(clippy::too_many_arguments)]
/// What decides how many alignments the seeding stage does.
///
/// `static` aligns everything the prefilter passed, so `max_seqs` bounds the
/// work directly. `prog` aligns from `prog_n` upward while the hit fraction
/// holds above `prog_f`. A `None` leaves nail on its own default.
pub struct Seeding<'a> {
    pub mmseqs_s: &'a str,
    pub mode: &'a str,
    /// `--mmseqs-max-seqs`: how much of the prefilter is allowed through.
    pub max_seqs: Option<usize>,
    /// `--prog-n`: the first round's alignments per query. prog mode only.
    pub prog_n: Option<usize>,
    /// `--prog-f`: the hit fraction that keeps it going. prog mode only.
    pub prog_f: Option<f64>,
}

impl<'a> Seeding<'a> {
    /// Seeding at a sensitivity and mode, every other knob left to nail.
    pub fn new(mmseqs_s: &'a str, mode: &'a str) -> Seeding<'a> {
        Seeding {
            mmseqs_s,
            mode,
            max_seqs: None,
            prog_n: None,
            prog_f: None,
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn seed(
    nail: &Path,
    mmseqs: &Path,
    query_hmm: &Path,
    target: &Path,
    shard_name: &str,
    seeds_out: &Path,
    dirs: &Dirs,
    threads: usize,
    seeding: &Seeding,
    fields: &[(&str, String)],
) -> Step {
    let mut cmd = Cmd::new(nail)
        .sub("search")
        .arg("--mmseqs-path", mmseqs)
        .arg("-t", threads)
        .arg("--tmp-dir", dirs.tmp.join("seeding"))
        .arg("--mmseqs-s", seeding.mmseqs_s)
        .arg("--seed-mode", seeding.mode);

    if let Some(n) = seeding.max_seqs {
        cmd = cmd.arg("--mmseqs-max-seqs", n);
    }
    if let Some(n) = seeding.prog_n {
        cmd = cmd.arg("--prog-n", n);
    }
    if let Some(f) = seeding.prog_f {
        cmd = cmd.arg("--prog-f", f);
    }

    let cmd = cmd
        .arg("--seeds-out", seeds_out)
        .flag("--only-seed")
        .flag("--allow-overwrite")
        .path(query_hmm)
        .path(target)
        .field(manifest::SHARD, shard_name);

    let cmd = fields
        .iter()
        .fold(cmd, |cmd, (key, value)| cmd.field(*key, value));

    Step::serial([cmd]).name("seeds").cores(threads)
}

/// One mmseqs search over one target, and the conversion that writes its table.
///
/// mmseqs takes every core it can find unless told otherwise, so the
/// conversion is held to the same count as the search rather than left to help
/// itself while something else is being timed.
pub struct Mmseqs<'a> {
    pub bin: &'a Path,
    pub query_db: &'a Path,
    pub target_db: &'a Path,
    /// Where the alignments land. Its parent has to exist already.
    pub aln_db: PathBuf,
    /// mmseqs' own scratch for the search.
    pub work: PathBuf,
    /// The hit table, in blast tabular form.
    pub out: PathBuf,
    pub threads: usize,
    /// `-s`, for a run that sweeps sensitivity. mmseqs' own default otherwise.
    pub s: Option<String>,
    /// `--max-seqs`. mmseqs' own default of 300 otherwise, which loses hits
    /// nail's seeding keeps.
    pub max_seqs: Option<usize>,
    /// `-e`. A benchmark drawing an ROC curve wants every hit reported.
    pub evalue: &'a str,
}

/// The two commands one mmseqs run takes.
///
/// Named rather than a pair, so a caller can field them separately: what the
/// run cost is the search, and the conversion is not part of it.
pub struct MmseqsCmds {
    pub search: Cmd,
    pub convert: Cmd,
}

impl Mmseqs<'_> {
    pub fn cmds(&self) -> MmseqsCmds {
        let mut search = Cmd::new(self.bin)
            .name("search")
            .sub("search")
            .arg("--threads", self.threads)
            .arg("-k", MMSEQS_K);

        if let Some(s) = &self.s {
            search = search.arg("-s", s);
        }
        if let Some(max) = self.max_seqs {
            search = search.arg("--max-seqs", max);
        }

        MmseqsCmds {
            search: search
                .arg("-e", self.evalue)
                .path(self.query_db)
                .path(self.target_db)
                .path(&self.aln_db)
                .path(&self.work),
            convert: Cmd::new(self.bin)
                .name("convertalis")
                .sub("convertalis")
                .arg("--threads", self.threads)
                .arg("--format-mode", 0)
                .path(self.query_db)
                .path(self.target_db)
                .path(&self.aln_db)
                .path(&self.out)
                .field(manifest::STAGE, CONVERT),
        }
    }
}

/// How many ways a query splits for hmmer, given a thread budget.
pub fn jobs(threads: usize) -> usize {
    (threads / HMMER_CPU).max(1)
}

/// The three binaries every pipeline needs, checked before anything runs.
pub struct Bins {
    pub nail: PathBuf,
    pub mmseqs: PathBuf,
    pub hmmsearch: PathBuf,
}

impl Bins {
    pub fn find() -> anyhow::Result<Bins> {
        Ok(Bins {
            nail: crate::tools::nail().context("nail")?,
            mmseqs: crate::tools::mmseqs().context("mmseqs")?,
            hmmsearch: crate::tools::hmmsearch().context("hmmsearch")?,
        })
    }
}
