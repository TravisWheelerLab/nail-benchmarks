//! Searches the family profiles against the originals of the decoys any run
//! has to settle, with hmmsearch as the one judge for every tool.
//!
//! A decoy is a reversed TrEMBL sequence, and the reversal of a true homolog
//! scores well against the homolog's family. So a (query, decoy) pair whose
//! family also hits the decoy's original is a reject rather than a false
//! positive. Only the decoys a run ranked at or above its worst true pair can
//! move its ROC, so the originals searched are the union of those over every
//! run.

use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, ensure};
use clap::Parser;

use libsail::collection::{Indexable, Iterable};
use libsail::seq::fasta::{DEFAULT_LINE_WIDTH, IndexedFasta};
use michi::{Closure, Cmd, PipelineBuilder, Progress, Step, Table};
use util::ledger;
use util::split::Kind;

use util::search::{Dirs, HMMER_CPU, Split, jobs, tag};

use crate::parse::{self, Benchmark};
use crate::run::{EVALUE, MODE, PRF};

/// The judge's run name, and the tool it records.
pub const JUDGE: &str = "judge";
const JUDGE_TOOL: &str = "hmmer";

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
        args.threads.is_multiple_of(HMMER_CPU),
        "--threads needs to be a multiple of {HMMER_CPU} (for hmmer)"
    );

    let hmmsearch = util::tools::hmmsearch()?;

    let inp = crate::Inputs::open(&paths.set)?;
    let bm = Benchmark::new(&inp.truth)?;

    // (tool, mode) -> the decoys its runs have to settle, kept apart only so
    // tested.tbl can say what each tool contributed to the union
    let mut by_tool: BTreeMap<(String, String), HashSet<String>> = BTreeMap::new();

    for run in parse::runs(&paths.search())? {
        let hits = run.hits()?;
        by_tool
            .entry((run.tool.clone(), run.mode.clone()))
            .or_default()
            .extend(
                parse::tested(&hits, &bm, &run.mode)
                    .into_iter()
                    .map(str::to_string),
            );
    }

    ensure!(
        !by_tool.is_empty(),
        "no finished runs in {}",
        paths.search().display()
    );

    let decoys: HashSet<String> = by_tool.values().flatten().cloned().collect();
    for ((tool, mode), settle) in &by_tool {
        println!("{tool}.{mode}: {} decoys to settle", settle.len());
    }
    println!("{} decoys to settle in all", decoys.len());

    // the judge's E-values are on the scale of the search they settle rather
    // than of the smaller file it reads, so a reject means the same thing
    // however many decoys the union holds
    let z = IndexedFasta::open(&inp.target_fa)
        .with_context(|| format!("failed to open {}", inp.target_fa.display()))?
        .len();

    let mut dirs = Dirs::new(paths.reject(), paths.tmp.join("reject"));
    if let Some(tmp) = args.tmp {
        dirs.tmp = tmp;
    }

    let target = dirs.tmp.join("originals.fa");
    let parts = dirs.tmp.join("parts");

    let subset = {
        let (originals, target) = (inp.originals.clone(), target.clone());
        let tested = dirs.root.join(parse::TESTED);
        let decoys = decoys.clone();

        let mut rows: Vec<Vec<String>> = by_tool
            .iter()
            .map(|((tool, mode), settle)| {
                vec![tool.clone(), mode.clone(), settle.len().to_string()]
            })
            .collect();
        rows.push(vec![
            "all".to_string(),
            "-".to_string(),
            decoys.len().to_string(),
        ]);

        move || -> anyhow::Result<()> {
            write_originals(&originals, &decoys, &target)?;
            let style = tabl::Style::default()
                .marker(tabl::Marker::Indent)
                .trailing(tabl::Trailing::Keep);
            let mut table =
                tabl::Table::new(tabl::Schema::new(["tool", "mode", "decoys"]).style(style));
            for row in rows {
                table.row(row);
            }
            table
                .write(&tested)
                .with_context(|| format!("failed to write {}", tested.display()))
        }
    };

    let mut pl = PipelineBuilder::new()
        .step(dirs.clean())
        .step(Cmd::new("mkdir").name("dirs").flag("-p").path(&dirs.tmp))
        .step(Step::from_closures([Closure::new("subset", subset)]).name("subset"));

    if !decoys.is_empty() {
        let split = Split::new(&inp.query_hmm, Kind::Hmm, &parts, jobs(args.threads));
        let fields = [(MODE, PRF.to_string())];

        let search = split.parts().iter().enumerate().map(|(i, part)| {
            let cmd = Cmd::new(&hmmsearch)
                .name(i.to_string())
                .arg("--cpu", HMMER_CPU)
                .arg("-Z", z)
                .arg("--tblout", parts.join(format!("{i}.tbl")))
                .arg("-E", EVALUE)
                .path(part)
                .path(&target);

            tag(cmd, JUDGE, JUDGE_TOOL, &fields)
        });

        let cat = (0..split.parts().len())
            .fold(Cmd::new("cat").name("tbl"), |cmd, i| {
                cmd.path(parts.join(format!("{i}.tbl")))
            })
            .stdout_to(dirs.table(JUDGE, ""));

        pl = pl
            .step(split.step(JUDGE, &[]))
            .step(
                Step::batched(split.parts().len(), search)
                    .name(JUDGE)
                    .cores(HMMER_CPU),
            )
            .step(
                Step::serial([tag(cat, JUDGE, JUDGE_TOOL, &fields)]).name(format!("cat.{JUDGE}")),
            );
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

/// Every original named in `decoys`, written into `to`.
fn write_originals(originals: &Path, decoys: &HashSet<String>, to: &Path) -> anyhow::Result<()> {
    let fa = IndexedFasta::open(originals)
        .with_context(|| format!("failed to open {}", originals.display()))?;

    let file = File::create(to).with_context(|| format!("failed to create {}", to.display()))?;
    let mut writer = BufWriter::new(file);

    for rec in fa.iter() {
        if decoys.contains(rec.name_str()?) {
            rec.write_to(&mut writer, DEFAULT_LINE_WIDTH)?;
        }
    }

    writer.flush()?;
    Ok(())
}
