//! What each search tool has to be told to do.
//!
//! Nothing here runs a process or touches the filesystem, apart from HMMER,
//! which has to read its query set to know how many ways to split it. Tools
//! describe work; [`crate::measure`] and [`crate::batch`] carry it out.

use std::path::{Path, PathBuf};

use anyhow::{bail, Context};

use crate::config::Run;
use crate::exec::Cmd;
use crate::{Paths, Search};
use bioio::split::{self, Kind};

/// The tools a benchmark can invoke.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Tool {
    Nail,
    Hmmer,
    Mmseqs,
    Blast,
    Last,
    Diamond,
}

/// A database or index built once and reused by every search that needs it.
pub struct Setup {
    /// The file whose existence proves this step already ran.
    pub marker: PathBuf,
    /// Directories to create before the commands run.
    pub dirs: Vec<PathBuf>,
    pub cmds: Vec<Cmd>,
    /// Named in the error when a command fails.
    pub what: String,
}

/// How the commands that make up a measurement are run.
pub enum Shape {
    One(Cmd),
    /// All at once, timings summed. HMMER scales poorly past a few threads, so
    /// its query set is split.
    Together(Vec<Cmd>),
    /// In turn, timings summed. psiblast takes one alignment at a time.
    Each(Vec<Cmd>),
}

/// Bookkeeping once the measured commands have finished.
pub enum After {
    /// Join split outputs into the one table the analysis expects.
    Concat { parts: Vec<PathBuf>, into: PathBuf },
    /// Move a file the tool wrote where it wanted it. Missing sources are
    /// ignored, since a tool may not produce one on every run.
    Move { from: PathBuf, to: PathBuf },
    /// A command that is not part of the measurement.
    Run { cmd: Cmd, what: String },
    Remove(PathBuf),
}

/// Everything one (run, search) pair does.
pub struct Work {
    /// Directories to create before the search runs, for tools that will not
    /// make their own.
    pub dirs: Vec<PathBuf>,
    pub search: Shape,
    pub after: Vec<After>,
}

impl Tool {
    pub fn parse(name: &str) -> anyhow::Result<Tool> {
        Ok(match name {
            "nail" => Tool::Nail,
            "hmmer" => Tool::Hmmer,
            "mmseqs" => Tool::Mmseqs,
            "blast" => Tool::Blast,
            "last" => Tool::Last,
            "diamond" => Tool::Diamond,
            other => bail!(
                "unknown tool {other:?}; known tools: nail, hmmer, mmseqs, blast, last, diamond"
            ),
        })
    }

    pub fn name(self) -> &'static str {
        match self {
            Tool::Nail => "nail",
            Tool::Hmmer => "hmmer",
            Tool::Mmseqs => "mmseqs",
            Tool::Blast => "blast",
            Tool::Last => "last",
            Tool::Diamond => "diamond",
        }
    }

    /// Whether this tool wants a working directory of its own, emptied before
    /// each search.
    pub fn uses_scratch(self) -> bool {
        matches!(self, Tool::Nail | Tool::Hmmer | Tool::Mmseqs)
    }

    /// Databases this tool needs before it can search. Nail and HMMER read
    /// their inputs directly and need none.
    pub fn setup(self, search: &Search, paths: &Paths) -> anyhow::Result<Vec<Setup>> {
        match self {
            Tool::Nail | Tool::Hmmer => Ok(Vec::new()),
            Tool::Mmseqs => mmseqs_setup(search, paths),
            Tool::Blast => Ok(vec![makedb(
                paths.tool("makeblastdb")?,
                blast_db(search, paths),
                "pdb",
                |bin, db| {
                    Cmd::new(bin)
                        .arg("-in")
                        .path(&search.target)
                        .arg("-dbtype")
                        .arg("prot")
                        .arg("-out")
                        .path(db)
                },
                "makeblastdb",
            )]),
            Tool::Last => Ok(vec![makedb(
                paths.tool("lastdb")?,
                last_db(search, paths),
                "prj",
                |bin, db| Cmd::new(bin).arg("-p").path(db).path(&search.target),
                "lastdb",
            )]),
            Tool::Diamond => Ok(vec![makedb(
                paths.tool("diamond")?,
                diamond_db(search, paths),
                "dmnd",
                |bin, db| {
                    Cmd::new(bin)
                        .arg("makedb")
                        .arg("--in")
                        .path(&search.target)
                        .arg("--db")
                        .path(db)
                },
                "diamond makedb",
            )]),
        }
    }

    /// The commands one run of one search consists of.
    pub fn work(self, run: &Run, search: &Search, paths: &Paths) -> anyhow::Result<Work> {
        match self {
            Tool::Nail => nail(run, search, paths),
            Tool::Hmmer => hmmer(run, search, paths),
            Tool::Mmseqs => mmseqs(run, search, paths),
            Tool::Blast => blast(run, search, paths),
            Tool::Last => last(run, search, paths),
            Tool::Diamond => diamond(run, search, paths),
        }
    }
}

// ------------------------------------------------------------------- shared

/// Which side of a tool to invoke. Every tool here has a profile mode and a
/// sequence mode; each maps `prf`/`seq` onto its own inputs and binaries.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Query {
    Profile,
    Sequence,
}

fn query_kind(run: &Run) -> anyhow::Result<Query> {
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

fn sequence_only(run: &Run, tool: &str) -> anyhow::Result<()> {
    if query_kind(run)? != Query::Sequence {
        bail!(
            "{tool} has no profile mode; use query = \"seq\" in run {:?}",
            run.name
        );
    }
    Ok(())
}

/// A query artifact the run needs, or a message naming what the benchmark did
/// not provide.
fn need<'a>(path: &'a Option<PathBuf>, what: &str) -> anyhow::Result<&'a Path> {
    path.as_deref().with_context(|| {
        format!("this benchmark provides no {what}, which the requested tool and query mode need")
    })
}

/// A stable directory name for a path, so databases derived from a shared query
/// set are built once rather than once per shard.
fn slug(path: &Path) -> String {
    path.to_string_lossy()
        .chars()
        .map(|c| if c.is_alphanumeric() { c } else { '_' })
        .collect::<String>()
        .trim_matches('_')
        .to_string()
}

/// The single-command shape the three sequence-only tools share.
fn makedb(
    bin: PathBuf,
    db: PathBuf,
    marker_ext: &str,
    build: impl FnOnce(PathBuf, &Path) -> Cmd,
    what: &str,
) -> Setup {
    Setup {
        marker: db.with_extension(marker_ext),
        dirs: db.parent().map(Path::to_path_buf).into_iter().collect(),
        cmds: vec![build(bin, &db)],
        what: what.to_string(),
    }
}

// --------------------------------------------------------------------- nail

fn nail(run: &Run, search: &Search, paths: &Paths) -> anyhow::Result<Work> {
    let query = match query_kind(run)? {
        Query::Profile => need(&search.hmm, "an hmm profile file")?,
        Query::Sequence => need(&search.fasta, "a query fasta file")?,
    };

    let scratch = paths.scratch(Tool::Nail, run, search);

    let cmd = Cmd::new(paths.tool("nail")?)
        .arg("search")
        .arg("--mmseqs-path")
        .path(paths.tool("mmseqs")?)
        .arg("-t")
        .arg(run.threads.to_string())
        .arg("--tmp-dir")
        .path(&scratch)
        .arg("--tbl-out")
        .path(paths.out(run, search, "tbl"))
        .args(run.args.clone())
        .path(query)
        .path(&search.target)
        .stderr_to(paths.log(run, search));

    Ok(Work {
        dirs: Vec::new(),
        search: Shape::One(cmd),
        // nail drops seeds in its tmp dir; keep them beside the other outputs
        after: vec![After::Move {
            from: scratch.join("seeds.tsv"),
            to: paths.out(run, search, "seeds"),
        }],
    })
}

// -------------------------------------------------------------------- hmmer

fn hmmer(run: &Run, search: &Search, paths: &Paths) -> anyhow::Result<Work> {
    let (program, query, split_kind) = match query_kind(run)? {
        Query::Profile => ("hmmsearch", need(&search.hmm, "an hmm profile file")?, Kind::Hmm),
        Query::Sequence => ("phmmer", need(&search.fasta, "a query fasta file")?, Kind::Fasta),
    };

    let program = paths.tool(program)?;
    let tbl = paths.out(run, search, "tbl");
    // per-domain output as well as per-sequence, since analysis needs domain
    // scores
    let dom = paths.out(run, search, "domtbl");
    let log = paths.log(run, search);

    let invoke = |cpus: usize, query: &Path, tbl: &Path, dom: &Path| {
        Cmd::new(&program)
            .arg("--cpu")
            .arg(cpus.to_string())
            .args(run.args.clone())
            .arg("-o")
            .arg("/dev/null")
            .arg("--tblout")
            .path(tbl)
            .arg("--domtblout")
            .path(dom)
            .path(query)
            .path(&search.target)
            .stderr_to(&log)
    };

    // hmmsearch/phmmer scale poorly past a few threads, so the query set is
    // split and run as several concurrent processes when threads_per is set
    let splits = match run.threads_per {
        Some(per) if per < run.threads => run.threads / per,
        _ => 1,
    };

    if splits <= 1 {
        return Ok(Work {
            dirs: Vec::new(),
            search: Shape::One(invoke(run.threads, query, &tbl, &dom)),
            after: Vec::new(),
        });
    }

    // the one place a tool reads anything: how many parts a query set actually
    // yields is not known until it has been indexed, and the count decides how
    // many processes to spawn
    let scratch = paths.scratch(Tool::Hmmer, run, search);
    let parts = split::write_splits(query, split_kind, splits, scratch.join("query"))?;
    let per = run.threads_per.expect("splits > 1 implies threads_per");

    let mut cmds = Vec::with_capacity(parts.len());
    let mut part_tbls = Vec::with_capacity(parts.len());
    let mut part_doms = Vec::with_capacity(parts.len());

    for part in &parts {
        let part_tbl = part.with_extension("tbl");
        let part_dom = part.with_extension("domtbl");
        cmds.push(invoke(per, part, &part_tbl, &part_dom));
        part_tbls.push(part_tbl);
        part_doms.push(part_dom);
    }

    Ok(Work {
        dirs: Vec::new(),
        search: Shape::Together(cmds),
        after: vec![
            After::Concat { parts: part_tbls, into: tbl },
            After::Concat { parts: part_doms, into: dom },
        ],
    })
}

// ------------------------------------------------------------------- mmseqs

fn mmseqs_target_db(search: &Search, paths: &Paths) -> PathBuf {
    paths
        .tmp
        .join(format!("mmseqs/target-{}/db", slug(&search.target)))
}

fn mmseqs_query_db(search: &Search, paths: &Paths, query: Query) -> anyhow::Result<PathBuf> {
    let (kind, source) = match query {
        Query::Profile => ("prf", need(&search.sto, "a stockholm alignment file")?),
        Query::Sequence => ("seq", need(&search.fasta, "a query fasta file")?),
    };
    Ok(paths.tmp.join(format!("mmseqs/query-{kind}-{}/db", slug(source))))
}

fn mmseqs_setup(search: &Search, paths: &Paths) -> anyhow::Result<Vec<Setup>> {
    let mmseqs = paths.tool("mmseqs")?;
    let mut out = Vec::new();

    // one target db per shard; query dbs are keyed by source path, so a query
    // set shared across shards is converted once rather than once per shard
    let tdb = mmseqs_target_db(search, paths);
    out.push(Setup {
        dirs: tdb.parent().map(Path::to_path_buf).into_iter().collect(),
        cmds: vec![Cmd::new(&mmseqs)
            .arg("createdb")
            .path(&search.target)
            .path(&tdb)],
        marker: tdb,
        what: "mmseqs createdb (target)".to_string(),
    });

    if search.fasta.is_some() {
        let qdb = mmseqs_query_db(search, paths, Query::Sequence)?;
        out.push(Setup {
            dirs: qdb.parent().map(Path::to_path_buf).into_iter().collect(),
            cmds: vec![Cmd::new(&mmseqs)
                .arg("createdb")
                .path(need(&search.fasta, "a query fasta file")?)
                .path(&qdb)],
            marker: qdb,
            what: "mmseqs createdb (query)".to_string(),
        });
    }

    if let Some(sto) = &search.sto {
        let qdb = mmseqs_query_db(search, paths, Query::Profile)?;
        let dir = qdb.parent().expect("db path has a parent").to_path_buf();
        let msa_db = dir.join("msaDB");

        out.push(Setup {
            cmds: vec![
                Cmd::new(&mmseqs)
                    .arg("convertmsa")
                    .path(sto)
                    .path(&msa_db)
                    .arg("--identifier-field")
                    .arg("0"),
                Cmd::new(&mmseqs)
                    .arg("msa2profile")
                    .path(&msa_db)
                    .path(&qdb)
                    .arg("--match-mode")
                    .arg("1"),
            ],
            dirs: vec![dir],
            marker: qdb,
            what: "mmseqs msa2profile".to_string(),
        });
    }

    Ok(out)
}

fn mmseqs(run: &Run, search: &Search, paths: &Paths) -> anyhow::Result<Work> {
    let mmseqs = paths.tool("mmseqs")?;
    let qdb = mmseqs_query_db(search, paths, query_kind(run)?)?;
    let tdb = mmseqs_target_db(search, paths);
    let log = paths.log(run, search);

    // mmseqs aborts rather than overwrite an existing alignment db, so both it
    // and the working directory live in scratch, which is emptied per search
    let scratch = paths.scratch(Tool::Mmseqs, run, search);
    let adb = scratch.join("align/db");
    let work = scratch.join("work");

    let cmd = Cmd::new(&mmseqs)
        // mmseqs reports failures on stdout, not stderr, so it shares the run
        // log; the log is discarded when the run succeeds
        .stdout_append(&log)
        .arg("search")
        .path(&qdb)
        .path(&tdb)
        .path(&adb)
        .path(&work)
        .arg("--threads")
        .arg(run.threads.to_string())
        .args(run.args.clone())
        .stderr_to(&log);

    let convert = Cmd::new(&mmseqs)
        .stdout_append(&log)
        .arg("convertalis")
        .path(&qdb)
        .path(&tdb)
        .path(&adb)
        .path(paths.out(run, search, "tbl"))
        .arg("--format-mode")
        .arg("0")
        .stderr_to(&log);

    Ok(Work {
        // mmseqs writes its alignment db and its working files into
        // directories it expects to find already there
        dirs: vec![
            adb.parent().expect("db path has a parent").to_path_buf(),
            work.clone(),
        ],
        search: Shape::One(cmd),
        after: vec![
            After::Run { cmd: convert, what: "mmseqs convertalis".to_string() },
            After::Remove(work),
        ],
    })
}

// -------------------------------------------------------------------- blast

fn blast_db(search: &Search, paths: &Paths) -> PathBuf {
    paths
        .tmp
        .join(format!("blast/{}/target_db", slug(&search.target)))
}

fn blast(run: &Run, search: &Search, paths: &Paths) -> anyhow::Result<Work> {
    let tbl = paths.out(run, search, "tbl");
    let log = paths.log(run, search);
    let db = blast_db(search, paths);

    match query_kind(run)? {
        Query::Sequence => Ok(Work {
            dirs: Vec::new(),
            search: Shape::One(
                Cmd::new(paths.tool("blastp")?)
                    .arg("-query")
                    .path(need(&search.fasta, "a query fasta file")?)
                    .arg("-db")
                    .path(&db)
                    .arg("-out")
                    .path(&tbl)
                    .arg("-outfmt")
                    .arg("6")
                    .arg("-num_threads")
                    .arg(run.threads.to_string())
                    .args(run.args.clone())
                    .stderr_to(&log),
            ),
            after: Vec::new(),
        }),

        // psiblast takes one alignment at a time, so a profile run is one
        // invocation per family with output collected into a single table
        Query::Profile => {
            let psiblast = paths.tool("psiblast")?;
            let afa_dir = need(&search.afa, "a directory of aligned fasta")?;

            let mut msas: Vec<PathBuf> = std::fs::read_dir(afa_dir)
                .with_context(|| format!("failed to read {}", afa_dir.display()))?
                .filter_map(|e| e.ok().map(|e| e.path()))
                .filter(|p| p.extension().is_some_and(|e| e == "afa"))
                .collect();
            msas.sort();

            if msas.is_empty() {
                bail!("no .afa files in {}", afa_dir.display());
            }

            let cmds = msas
                .iter()
                .enumerate()
                .map(|(i, msa)| {
                    let cmd = Cmd::new(&psiblast)
                        .arg("-in_msa")
                        .path(msa)
                        .arg("-db")
                        .path(&db)
                        .arg("-outfmt")
                        .arg("6")
                        .arg("-num_threads")
                        .arg(run.threads.to_string())
                        .arg("-comp_based_stats")
                        .arg("1")
                        .arg("-num_iterations")
                        .arg("1")
                        .args(run.args.clone())
                        .stderr_to(&log);

                    // the first invocation truncates whatever an earlier run
                    // left behind; the rest append to it
                    if i == 0 {
                        cmd.stdout_to(&tbl)
                    } else {
                        cmd.stdout_append(&tbl)
                    }
                })
                .collect();

            Ok(Work {
                dirs: Vec::new(),
                search: Shape::Each(cmds),
                after: Vec::new(),
            })
        }
    }
}

// --------------------------------------------------------------------- last

fn last_db(search: &Search, paths: &Paths) -> PathBuf {
    paths
        .tmp
        .join(format!("last/{}/target_db", slug(&search.target)))
}

fn last(run: &Run, search: &Search, paths: &Paths) -> anyhow::Result<Work> {
    sequence_only(run, "last")?;

    // lastal writes its table to stdout
    let cmd = Cmd::new(paths.tool("lastal")?)
        .path(last_db(search, paths))
        .path(need(&search.fasta, "a query fasta file")?)
        .arg("-f")
        .arg("BlastTab")
        .arg("-P")
        .arg(run.threads.to_string())
        .args(run.args.clone())
        .stdout_to(paths.out(run, search, "tbl"))
        .stderr_to(paths.log(run, search));

    Ok(Work {
        dirs: Vec::new(),
        search: Shape::One(cmd),
        after: Vec::new(),
    })
}

// ------------------------------------------------------------------ diamond

fn diamond_db(search: &Search, paths: &Paths) -> PathBuf {
    paths
        .tmp
        .join(format!("diamond/{}/target_db", slug(&search.target)))
}

fn diamond(run: &Run, search: &Search, paths: &Paths) -> anyhow::Result<Work> {
    sequence_only(run, "diamond")?;

    let cmd = Cmd::new(paths.tool("diamond")?)
        .arg("blastp")
        .arg("--query")
        .path(need(&search.fasta, "a query fasta file")?)
        .arg("--db")
        .path(diamond_db(search, paths))
        .arg("--out")
        .path(paths.out(run, search, "tbl"))
        .arg("--outfmt")
        .arg("6")
        .arg("--threads")
        .arg(run.threads.to_string())
        .args(run.args.clone())
        .stderr_to(paths.log(run, search));

    Ok(Work {
        dirs: Vec::new(),
        search: Shape::One(cmd),
        after: Vec::new(),
    })
}
