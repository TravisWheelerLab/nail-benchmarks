//! The steps the run is assembled out of.
//!
//! These return [`Step`]s and [`Cmd`]s for the pipeline to compose. Nothing
//! here owns a pipeline or decides what a run measures.
//!
//! Every search command carries the fields `parse` reads back: `name` for the
//! column it becomes and the stem of its results file, `tool` for how to read
//! that file, and `mode` for which query it searched. There is one target set,
//! so nothing carries a shard.

use std::path::{Path, PathBuf};

use anyhow::Context;

use bench::manifest;
use bioio::split::{self, Kind};
use pail::{Closure, Cmd, Step};

use crate::inputs::Inputs;

/// hmmsearch and phmmer scale poorly past a couple of threads, so the query set
/// is split threads/HMMER_CPU ways and the parts run at the same time.
pub const HMMER_CPU: usize = 2;

/// Every tool but blast reports down to here, so they are comparable. blast is
/// deliberately left at its default: matching this makes it dramatically slower
/// for no extra recall.
pub const EVALUE: &str = "1e9";
pub const MAX_SEQS: &str = "2000";

/// Which query a run searched, which is the axis this benchmark exists for.
pub const MODE: &str = "mode";
pub const PRF: &str = "prf";
pub const SEQ: &str = "seq";

/// Where one run's output lives.
pub struct Dirs {
    pub root: PathBuf,
    pub results: PathBuf,
    pub tmp: PathBuf,
}

impl Dirs {
    pub fn new(set: &Inputs) -> Dirs {
        let root = set.run_dir();
        Dirs {
            results: root.join("results"),
            tmp: root.join("tmp"),
            root,
        }
    }

    /// The hit table a run writes. No shard: there is one target set, so a run
    /// is named by what searched it and nothing else.
    pub fn table(&self, name: &str) -> PathBuf {
        manifest::table_path(&self.results, name, "")
    }

    /// Everything a run leaves behind, cleared at the front of the pipeline
    /// rather than before it.
    ///
    /// A stale table from an earlier run reads as a run that simply found
    /// less, and mmseqs refuses to overwrite an existing alignment db. Doing it
    /// as a step is what keeps `--dry-run` from touching the disk.
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
}

/// A query set cut into parts for hmmer to search in parallel.
///
/// The parts are named before they exist: `write_splits` files them by index,
/// so the batch can be written without waiting to see what the split produced.
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
        let ext = match kind {
            Kind::Hmm => "hmm",
            Kind::Fasta => "fa",
        };

        Split {
            query: query.into(),
            kind,
            parts: (0..jobs).map(|i| dir.join(format!("{i}.{ext}"))).collect(),
            dir,
        }
    }

    /// Rust in place of a command, so it is a closure step rather than
    /// something that happens while the pipeline is being built. Whatever a
    /// previous run left in there would be searched as if it belonged.
    pub fn step(&self, name: &str) -> Step {
        let (query, kind, dir, jobs) = (
            self.query.clone(),
            self.kind,
            self.dir.clone(),
            self.parts.len(),
        );

        Step::from_closures([Closure::new("split", move || {
            std::fs::remove_dir_all(&dir).ok();
            let written = split::write_splits(&query, kind, jobs, &dir)?;

            // empty bins are skipped, so a query with fewer records than parts
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
        .name(format!("{name}.split"))
    }
}

/// One hmmer-family run: the query's parts searched together, then their hits
/// gathered into one table.
pub struct Hmmer {
    pub search: Step,
    pub cat: Step,
}

/// What a run is called and what produced it: the three fields every search
/// command carries, passed together because they always travel together.
#[derive(Clone, Copy)]
pub struct Run<'a> {
    pub name: &'a str,
    pub tool: &'a str,
    pub mode: &'a str,
}

/// hmmsearch over profiles or phmmer over sequences -- the same shape either
/// way, which is why they share a builder.
pub fn hmmer(
    program: &Path,
    split: &Split,
    dirs: &Dirs,
    run: Run<'_>,
    target: &Path,
    threads: usize,
) -> Hmmer {
    let Run { name, tool, mode } = run;

    let parts = &split.parts;
    let scratch = split.dir.clone();

    let fields = |cmd: Cmd| {
        cmd.field(manifest::NAME, name)
            .field(manifest::TOOL, tool)
            .field(MODE, mode)
    };

    Hmmer {
        search: Step::batched(
            parts.len(),
            parts.iter().enumerate().map(|(i, part)| {
                fields(
                    Cmd::new(program)
                        .name(i.to_string())
                        .arg("--cpu", HMMER_CPU)
                        .arg("-E", EVALUE)
                        .arg("-o", "/dev/null")
                        .arg("--tblout", scratch.join(format!("{i}.tbl")))
                        .arg("--domtblout", scratch.join(format!("{i}.domtbl")))
                        .path(part)
                        .path(target),
                )
            }),
        )
        .name(name)
        // per command, not per step, so this asks for HMMER_CPU x parts, which
        // is --threads again
        .cores(HMMER_CPU.min(threads)),
        cat: Step::serial([fields(
            cat(
                (0..parts.len()).map(|i| scratch.join(format!("{i}.tbl"))),
                dirs.table(name),
            )
            .name("tbl"),
        )])
        .name(format!("{name}.cat")),
    }
}

/// One mmseqs search and the conversion that writes its table.
///
/// mmseqs aborts rather than overwrite an existing alignment db, so the pair is
/// preceded by a clean of its own scratch.
pub struct Mmseqs<'a> {
    pub bin: &'a Path,
    pub query_db: &'a Path,
    pub target_db: &'a Path,
    pub scratch: PathBuf,
    pub out: PathBuf,
    pub threads: usize,
    pub s: String,
}

/// The scratch management and the run itself, kept apart because only the run
/// is the run: the fields go on the two mmseqs calls, so what the column cost
/// is the search plus the conversion and not the mkdir before them.
pub struct MmseqsCmds {
    pub prep: Vec<Cmd>,
    pub search: [Cmd; 2],
}

impl Mmseqs<'_> {
    pub fn cmds(&self) -> MmseqsCmds {
        let aln_db = self.scratch.join("alnDB");
        let work = self.scratch.join("work");

        let prep = vec![
            Cmd::new("rm").name("clean").flag("-rf").path(&self.scratch),
            Cmd::new("mkdir")
                .name("dirs")
                .flag("-p")
                .path(aln_db.parent().expect("alnDB has a parent"))
                .path(&work),
        ];

        let search = [
            Cmd::new(self.bin)
                .name("search")
                .sub("search")
                .path(self.query_db)
                .path(self.target_db)
                .path(&aln_db)
                .path(&work)
                .arg("--threads", self.threads)
                .arg("-s", &self.s)
                .arg("--max-seqs", MAX_SEQS)
                .arg("-e", EVALUE),
            Cmd::new(self.bin)
                .name("convertalis")
                .sub("convertalis")
                .path(self.query_db)
                .path(self.target_db)
                .path(&aln_db)
                .path(&self.out)
                .arg("--format-mode", 0),
        ];

        MmseqsCmds { prep, search }
    }
}

/// There's no shell to expand a glob, so the parts get named one by one.
pub fn cat(parts: impl IntoIterator<Item = PathBuf>, into: PathBuf) -> Cmd {
    parts
        .into_iter()
        .fold(Cmd::new("cat"), |cmd, part| cmd.path(part))
        .stdout_to(into)
}

/// How many ways a query splits for hmmer, given a thread budget.
pub fn jobs(threads: usize) -> usize {
    (threads / HMMER_CPU).max(1)
}

/// Every binary the run needs, checked before anything starts.
///
/// All of them, up front: a sweep that dies two tools in has spent its wall
/// time on numbers that are no longer comparable to the ones it did not reach.
pub struct Bins {
    pub nail: PathBuf,
    pub mmseqs: PathBuf,
    pub hmmsearch: PathBuf,
    pub phmmer: PathBuf,
    pub blastp: PathBuf,
    pub psiblast: PathBuf,
    pub makeblastdb: PathBuf,
    pub lastal: PathBuf,
    pub lastdb: PathBuf,
    pub diamond: PathBuf,
}

impl Bins {
    pub fn find() -> anyhow::Result<Bins> {
        Ok(Bins {
            nail: tools::nail().context("nail")?,
            mmseqs: tools::mmseqs().context("mmseqs")?,
            hmmsearch: tools::hmmsearch().context("hmmsearch")?,
            phmmer: tools::phmmer().context("phmmer")?,
            blastp: tools::blastp().context("blastp")?,
            psiblast: tools::psiblast().context("psiblast")?,
            makeblastdb: tools::makeblastdb().context("makeblastdb")?,
            lastal: tools::lastal().context("lastal")?,
            lastdb: tools::lastdb().context("lastdb")?,
            diamond: tools::diamond().context("diamond")?,
        })
    }
}
