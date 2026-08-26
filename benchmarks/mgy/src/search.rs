//! The steps every pipeline is assembled out of.
//!
//! These return [`Step`]s and [`Cmd`]s for a pipeline to compose. Nothing here
//! owns a pipeline or decides what a run measures — a pipeline that wants
//! hmmer against a shard asks for those steps and puts them where it wants
//! them.
//!
//! Every search command carries the fields `parse` reads back: `name` for the
//! column it becomes, `tool` for how to read its table, and `shard` for which
//! target it ran against. Seeding carries `stage` instead, since it produces
//! pairs rather than scores and so is not a column.

use std::path::{Path, PathBuf};

use anyhow::Context;

use bioio::split::{self, Kind};
use pail::{Closure, Cmd, Step};

use bench::manifest;

/// hmmsearch doesn't scale past a couple of threads, so its query gets split
/// threads/HMMER_CPU ways and the parts run at the same time.
pub const HMMER_CPU: usize = 2;

/// Every tool reports down to here, so they can be compared.
pub const EVALUE: &str = "10";

/// Where one pipeline's output lives.
///
/// Everything a run produces -- hit tables, domain tables, the seed list --
/// lands in `results/`, told apart by name rather than by directory. `tmp/`
/// holds the scratch that gets thrown away.
pub struct Dirs {
    pub root: PathBuf,
    pub results: PathBuf,
    pub tmp: PathBuf,
}

impl Dirs {
    pub fn new(name: &str) -> Dirs {
        let root = crate::runs().join(name);
        Dirs {
            results: root.join("results"),
            tmp: root.join("tmp"),
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

    pub fn table(&self, name: &str, shard: &str) -> PathBuf {
        manifest::table_path(&self.results, name, shard)
    }

    pub fn seeds(&self, shard: &str) -> PathBuf {
        manifest::seeds_path(&self.results, shard)
    }
}

/// There's no shell to expand a glob, so the parts get named one by one.
pub fn cat(parts: impl IntoIterator<Item = PathBuf>, into: PathBuf) -> Cmd {
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
    dir: PathBuf,
    parts: Vec<PathBuf>,
}

impl Split {
    pub fn new(query: impl Into<PathBuf>, dir: impl Into<PathBuf>, jobs: usize) -> Split {
        let dir = dir.into();
        Split {
            query: query.into(),
            parts: (0..jobs).map(|i| dir.join(format!("{i}.hmm"))).collect(),
            dir,
        }
    }

    /// Rust in place of a command, so it is a closure step. Whatever a previous
    /// run left in there would be searched as if it belonged.
    pub fn step(&self) -> Step {
        let (query, dir, jobs) = (self.query.clone(), self.dir.clone(), self.parts.len());

        Step::from_closures([Closure::new("split", move || {
            std::fs::remove_dir_all(&dir).ok();
            let written = split::write_splits(&query, Kind::Hmm, jobs, &dir)?;

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
        })])
        .name("split")
    }
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
pub fn hmmer(
    hmmsearch: &Path,
    split: &Split,
    dirs: &Dirs,
    name: &str,
    shard_name: &str,
    target: &Path,
    extra: &[(&str, String)],
) -> Hmmer {
    let parts = &split.parts;
    let scratch = dirs.tmp.join("hmmer");

    let fields = |cmd: Cmd| {
        extra.iter().fold(
            cmd.field(manifest::NAME, name)
                .field(manifest::TOOL, "hmmer")
                .field(manifest::SHARD, shard_name),
            |cmd, (key, value)| cmd.field(*key, value),
        )
    };

    Hmmer {
        search: Step::batched(
            parts.len(),
            parts.iter().enumerate().map(|(i, part)| {
                fields(
                    Cmd::new(hmmsearch)
                        .name(i.to_string())
                        .arg("--cpu", HMMER_CPU)
                        .arg("--tblout", scratch.join(format!("{i}.tbl")))
                        .arg("--domtblout", scratch.join(format!("{i}.domtbl")))
                        .arg("-E", EVALUE)
                        .path(part)
                        .path(target),
                )
            }),
        )
        .name(format!("{name}.{shard_name}"))
        // per command, not per step, so this asks for HMMER_CPU x parts, which
        // is --threads again. a machine with a smaller pool than that won't
        // fail, it will just run fewer of the parts at once
        .cores(HMMER_CPU),
        cat: Step::serial([
            fields(
                cat(
                    (0..parts.len()).map(|i| scratch.join(format!("{i}.tbl"))),
                    manifest::table_path(&dirs.results, name, shard_name),
                )
                .name("tbl"),
            ),
            cat(
                (0..parts.len()).map(|i| scratch.join(format!("{i}.domtbl"))),
                manifest::dom_path(&dirs.results, name, shard_name),
            )
            .name("domtbl"),
        ])
        .name(format!("cat.{name}.{shard_name}")),
    }
}

/// One seeding pass, kept so later searches can replay it.
///
/// The seeds are per-pair, so `parse` reads them back as the `seeded` column —
/// which is what lets a pipeline ask where a hit was lost rather than only
/// whether it survived.
#[allow(clippy::too_many_arguments)]
pub fn seed(
    nail: &Path,
    mmseqs: &Path,
    query_hmm: &Path,
    target: &Path,
    shard_name: &str,
    dirs: &Dirs,
    threads: usize,
    mmseqs_s: &str,
    seed_mode: &str,
) -> Step {
    Step::serial([Cmd::new(nail)
        .sub("search")
        .arg("--mmseqs-path", mmseqs)
        .arg("-t", threads)
        .arg("--tmp-dir", dirs.tmp.join("seeding"))
        .arg("--mmseqs-s", mmseqs_s)
        .arg("--seed-mode", seed_mode)
        .arg("--seeds-out", dirs.seeds(shard_name))
        .flag("--only-seed")
        .flag("--allow-overwrite")
        .path(query_hmm)
        .path(target)
        .field(manifest::STAGE, "seed")
        .field(manifest::SHARD, shard_name)])
    .name("seeds")
    .cores(threads)
}

/// One mmseqs search over one target, and the conversion that writes its table.
///
/// Two commands rather than one, and both carry the run's fields: the search
/// does the work and the conversion writes the table, so what the column cost
/// is the two together.
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
}

impl Mmseqs<'_> {
    pub fn cmds(&self) -> [Cmd; 2] {
        let mut search = Cmd::new(self.bin)
            .name("search")
            .sub("search")
            .arg("--threads", self.threads);

        if let Some(s) = &self.s {
            search = search.arg("-s", s);
        }
        if let Some(max) = self.max_seqs {
            search = search.arg("--max-seqs", max);
        }

        [
            search
                .arg("-e", EVALUE)
                .path(self.query_db)
                .path(self.target_db)
                .path(&self.aln_db)
                .path(&self.work),
            Cmd::new(self.bin)
                .name("convertalis")
                .sub("convertalis")
                .arg("--threads", self.threads)
                .arg("--format-mode", 0)
                .path(self.query_db)
                .path(self.target_db)
                .path(&self.aln_db)
                .path(&self.out),
        ]
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
            nail: tools::nail().context("nail")?,
            mmseqs: tools::mmseqs().context("mmseqs")?,
            hmmsearch: tools::hmmsearch().context("hmmsearch")?,
        })
    }
}
