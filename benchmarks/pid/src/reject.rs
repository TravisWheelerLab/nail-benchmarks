//! Searches each tool against the originals of the decoys it has to settle.
//!
//! A decoy is a reversed TrEMBL sequence, and the reversal of a true homolog
//! scores well against the homolog's family. So a (query, decoy) pair whose
//! query also hits the decoy's original is a reject rather than a false
//! positive. Only the decoys a run ranked at or above its worst true pair can
//! move its ROC, so those are what a tool is searched against here: the union
//! over its runs, at the most sensitive setting any of them used.

use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, bail, ensure};
use clap::Parser;

use libsail::collection::Iterable;
use libsail::seq::fasta::{DEFAULT_LINE_WIDTH, IndexedFasta};
use michi::{Closure, Cmd, OnError, Output, PipelineBuilder, Progress, Step, Table};
use util::ledger;
use util::split::Kind;

use util::search::{Dirs, Mmseqs, Split, hmmer, jobs, nail, tag};

use crate::parse::{self, Benchmark};
use crate::run::{self, EVALUE, MODE, PRF, SEQ};

/// The preset diamond settles its decoys at, the more sensitive of the two
/// `run` sweeps.
const DIAMOND_PRESET: &str = "ultra-sensitive";

/// The prefilter cap nail and mmseqs settle their decoys at.
//
// wide open, as cutoffs' reject stage runs: at the search stage's cap an
// original ranked below a query's 2000th hit is never reported, and its decoy
// stays a false positive it was never tested for
const MAX_SEQS_OPEN: usize = 1_000_000_000;

/// (tool, mode) -> the highest -s its runs swept, and the decoys it has to
/// settle.
type Judges = BTreeMap<(String, String), (Option<f32>, HashSet<String>)>;

#[derive(Parser, Debug)]
pub struct Args {
    /// Which label of paths.toml to read. Omit to list them
    #[arg(long = "in", value_name = "label")]
    pub label: Option<String>,

    #[arg(short, long, default_value_t = 24)]
    pub threads: usize,

    #[arg(long)]
    pub tmp: Option<PathBuf>,

    #[arg(long)]
    pub dry_run: bool,
}

pub fn main(args: Args, paths: &crate::Paths) -> anyhow::Result<()> {
    ensure!(
        args.threads.is_multiple_of(util::search::HMMER_CPU),
        "--threads needs to be a multiple of {} (for hmmer)",
        util::search::HMMER_CPU
    );

    let nail_bin = util::tools::nail()?;
    let mmseqs = util::tools::mmseqs()?;
    let hmmsearch = util::tools::hmmsearch()?;
    let phmmer = util::tools::phmmer()?;
    let blastp_bin = util::tools::blastp()?;
    let psiblast_bin = util::tools::psiblast()?;
    let makeblastdb = util::tools::makeblastdb()?;
    let lastal_bin = util::tools::lastal()?;
    let lastdb = util::tools::lastdb()?;
    let diamond_bin = util::tools::diamond()?;

    let inp = crate::Inputs::open(&paths.set)?;
    let bm = Benchmark::new(&inp.truth)?;

    let mut judges: Judges = BTreeMap::new();

    for run in parse::runs(&paths.search())? {
        let hits = run.hits()?;
        let (s, decoys) = judges.entry(parse::judge(&run)).or_default();

        decoys.extend(
            parse::tested(&hits, &bm, &run.mode)
                .into_iter()
                .map(str::to_string),
        );

        if let Some(run_s) = run.params.get("s") {
            let run_s: f32 = run_s
                .parse()
                .with_context(|| format!("{} recorded -s {run_s:?}", run.name))?;
            *s = Some(s.map_or(run_s, |s| s.max(run_s)));
        }
    }

    ensure!(
        !judges.is_empty(),
        "no finished runs in {}",
        paths.search().display()
    );

    for ((tool, mode), (_, decoys)) in &judges {
        println!("{tool}.{mode}: {} decoys to settle", decoys.len());
    }

    let mut dirs = Dirs::new(paths.reject(), paths.tmp.join("reject"));
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let scratch = |tool: &str, mode: &str| dirs.tmp.join(format!("{tool}.{mode}"));
    let settling: Vec<_> = judges
        .iter()
        .filter(|(_, (_, decoys))| !decoys.is_empty())
        .collect();

    let subset = {
        let originals = inp.originals.clone();
        let tested = dirs.root.join(parse::TESTED);

        let wanted: Vec<(PathBuf, HashSet<String>)> = settling
            .iter()
            .map(|((tool, mode), (_, decoys))| {
                (scratch(tool, mode).join("originals.fa"), decoys.clone())
            })
            .collect();

        let rows: Vec<Vec<String>> = judges
            .iter()
            .map(|((tool, mode), (_, decoys))| {
                vec![tool.clone(), mode.clone(), decoys.len().to_string()]
            })
            .collect();

        move || -> anyhow::Result<()> {
            write_originals(&originals, &wanted)?;
            util::tbl::write(
                &tested,
                util::tbl::Table {
                    meta: "",
                    headers: &["tool", "mode", "decoys"].map(String::from),
                    rows: &rows,
                    ragged_last: false,
                },
            )
        }
    };

    let mut pl = PipelineBuilder::new()
        .step(dirs.clean())
        .step(settling.iter().fold(
            Cmd::new("mkdir").name("dirs").flag("-p"),
            |cmd, ((t, m), _)| cmd.path(scratch(t, m)),
        ));

    pl = pl.step(Step::from_closures([Closure::new("subset", subset)]).name("subset"));

    for ((tool, mode), (s, decoys)) in settling {
        let name = format!("{tool}.{mode}");
        let scratch = scratch(tool, mode);
        let target = scratch.join("originals.fa");
        let out = dirs.table(&name, "");
        let fields = [(MODE, mode.clone())];

        let query = match mode.as_str() {
            PRF => &inp.query_hmm,
            SEQ => &inp.query_fa,
            _ => bail!("{name}: no query for mode {mode:?}"),
        };

        let s = || s.with_context(|| format!("{name}: its runs recorded no -s"));

        match tool.as_str() {
            "nail" => {
                let s = format!("{:.1}", s()?);

                pl = pl.step(nail(
                    &nail_bin,
                    &mmseqs,
                    query,
                    &target,
                    &out,
                    &scratch.join("nail"),
                    args.threads,
                    EVALUE,
                    &[
                        ("--mmseqs-s", s.clone()),
                        ("--mmseqs-max-seqs", MAX_SEQS_OPEN.to_string()),
                    ],
                    &[],
                    &name,
                    &[(MODE, mode.clone()), ("s", s)],
                ));
            }

            "hmmer" | "phmmer" => {
                let (program, kind) = match tool.as_str() {
                    "hmmer" => (&hmmsearch, Kind::Hmm),
                    _ => (&phmmer, Kind::Fasta),
                };

                let split = Split::new(query, kind, scratch.join("parts"), jobs(args.threads));
                let searched = hmmer(
                    program, &split, &dirs, &name, tool, "", &target, EVALUE, false, &fields,
                );

                pl = pl
                    .step(split.step(&name, &[]))
                    .step(searched.search)
                    .step(searched.cat);
            }

            "mmseqs" => {
                let s = format!("{:.1}", s()?);
                let target_db = scratch.join("targetDB/targetDB");
                let query_db = scratch.join("queryDB/queryDB");
                let msa_db = scratch.join("msaDB/msaDB");

                let mut db = vec![
                    Cmd::new("mkdir")
                        .name("dirs")
                        .flag("-p")
                        .path(scratch.join("targetDB"))
                        .path(scratch.join("queryDB"))
                        .path(scratch.join("msaDB")),
                    Cmd::new(&mmseqs)
                        .name("createdb-target")
                        .sub("createdb")
                        .path(&target)
                        .path(&target_db),
                ];

                // the same query databases run builds, over the same files
                db.extend(match mode.as_str() {
                    PRF => vec![
                        Cmd::new(&mmseqs)
                            .name("convertmsa")
                            .sub("convertmsa")
                            .path(&inp.query_sto)
                            .path(&msa_db)
                            .arg("--identifier-field", 0),
                        Cmd::new(&mmseqs)
                            .name("msa2profile")
                            .sub("msa2profile")
                            .path(&msa_db)
                            .path(&query_db)
                            .arg("--match-mode", 1),
                    ],
                    _ => vec![
                        Cmd::new(&mmseqs)
                            .name("createdb-query")
                            .sub("createdb")
                            .path(query)
                            .path(&query_db),
                    ],
                });

                let cmds = Mmseqs {
                    bin: &mmseqs,
                    query_db: &query_db,
                    target_db: &target_db,
                    aln_db: scratch.join("alnDB"),
                    work: scratch.join("work"),
                    out,
                    threads: args.threads,
                    s: Some(s.clone()),
                    max_seqs: Some(MAX_SEQS_OPEN),
                    evalue: EVALUE,
                }
                .cmds();

                let fields = [(MODE, mode.clone()), ("s", s)];
                let searched =
                    [cmds.search, cmds.convert].map(|cmd| tag(cmd, &name, tool, &fields));

                pl = pl
                    .step(Step::serial(db).name(format!("{name}-db")))
                    .step(Step::serial(searched).name(name));
            }

            "blast" => {
                let db = scratch.join("blast/db");
                pl = pl.step(
                    Step::serial([
                        Cmd::new("mkdir")
                            .name("dirs")
                            .flag("-p")
                            .path(scratch.join("blast")),
                        Cmd::new(&makeblastdb)
                            .name("makeblastdb")
                            .arg("-in", &target)
                            .arg("-dbtype", "prot")
                            .arg("-out", &db),
                    ])
                    .name(format!("{name}-db")),
                );

                pl = match mode.as_str() {
                    SEQ => pl.step(run::blastp(
                        &blastp_bin,
                        query,
                        &db,
                        &out,
                        args.threads,
                        &name,
                        Some(decoys.len()),
                    )),
                    _ => {
                        let afa = run::alignments(&inp.afa)?;
                        pl.step(
                            Step::serial(afa.iter().enumerate().map(|(i, msa)| {
                                let cmd = run::psiblast(
                                    &psiblast_bin,
                                    msa,
                                    &db,
                                    args.threads,
                                    &name,
                                    Some(decoys.len()),
                                );
                                match i {
                                    0 => cmd.stdout_to(&out),
                                    _ => cmd.stdout(Output::Append(out.clone())),
                                }
                            }))
                            .name(name)
                            .on_error(OnError::Continue),
                        )
                    }
                };
            }

            "last" => {
                let db = scratch.join("last/db");
                pl = pl
                    .step(
                        Step::serial([
                            Cmd::new("mkdir")
                                .name("dirs")
                                .flag("-p")
                                .path(scratch.join("last")),
                            Cmd::new(&lastdb)
                                .name("lastdb")
                                .arg("-p", &db)
                                .path(&target),
                        ])
                        .name(format!("{name}-db")),
                    )
                    .step(run::lastal(
                        &lastal_bin,
                        &db,
                        query,
                        &out,
                        args.threads,
                        &name,
                    ));
            }

            "diamond" => {
                let db = scratch.join("diamond/db");
                pl = pl
                    .step(
                        Step::serial([
                            Cmd::new("mkdir")
                                .name("dirs")
                                .flag("-p")
                                .path(scratch.join("diamond")),
                            Cmd::new(&diamond_bin)
                                .name("diamond-makedb")
                                .sub("makedb")
                                .arg("--in", &target)
                                .arg("--db", &db),
                        ])
                        .name(format!("{name}-db")),
                    )
                    .step(run::diamond(
                        &diamond_bin,
                        query,
                        &db,
                        &out,
                        args.threads,
                        &name,
                        DIAMOND_PRESET,
                        Some(0),
                    ));
            }

            _ => bail!("{name}: no reject search for tool {tool:?}"),
        }
    }

    let pipeline = pl
        .stderr_dir(dirs.tmp.join("stderr"))
        .sink(Progress::new())
        .sink(Table::new(dirs.root.join("manifest.tbl")))
        .build()
        .context("failed to build the reject stage")?;

    if args.dry_run {
        pipeline.dry_run();
        return Ok(());
    }

    ledger::clear(&dirs.root);
    pipeline.run()?;
    ledger::record(&dirs.root)
}

/// Every original named in `wanted`, written into the file paired with it.
fn write_originals(originals: &Path, wanted: &[(PathBuf, HashSet<String>)]) -> anyhow::Result<()> {
    let fa = IndexedFasta::open(originals)
        .with_context(|| format!("failed to open {}", originals.display()))?;

    let mut writers = wanted
        .iter()
        .map(|(path, _)| {
            let file = File::create(path)
                .with_context(|| format!("failed to create {}", path.display()))?;
            Ok(BufWriter::new(file))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;

    for rec in fa.iter() {
        let name = rec.name_str()?;
        for ((_, names), writer) in wanted.iter().zip(&mut writers) {
            if names.contains(name) {
                rec.write_to(writer, DEFAULT_LINE_WIDTH)?;
            }
        }
    }

    for mut writer in writers {
        writer.flush()?;
    }

    Ok(())
}
