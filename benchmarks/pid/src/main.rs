//! Percent-identity benchmark: profmark-style ROC over Pfam families embedded
//! in a background of reversed TrEMBL sequences.
//!
//! The set is built by `build-set --in pid-toy|pid-real` and read from the
//! store like every other benchmark's. One label fills in every path, so a
//! subcommand never names a directory.

mod parse;
mod plot;
mod reject;
mod run;

use std::path::{Path, PathBuf};

use anyhow::Context;
use clap::{Parser, Subcommand};

use serde::Deserialize;
use util::set::{Set, shape};

#[derive(Parser)]
#[command(name = "pid", about = "percent-identity benchmark")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Search every tool against the benchmark.
    Run(run::Args),
    /// Search each tool against the originals of the decoys it ranked above
    /// its worst true pair.
    Reject(reject::Args),
    /// Turn results into the tables the plot scripts consume.
    #[command(subcommand)]
    Parse(parse::Cmd),
    /// Draw the figures from what parse wrote.
    Plot(plot::Args),
}

/// Where this benchmark reads and writes, out of `paths.toml`.
#[derive(Deserialize)]
pub struct Paths {
    /// The set to search.
    pub set: PathBuf,
    /// Where the results and the ledger go.
    pub run: PathBuf,
    /// Where the tables parse works out go.
    pub analysis: PathBuf,
    /// Where the figures go, which is outside the store: a pdf is read by a
    /// person rather than by another pipeline.
    pub figures: PathBuf,
    /// Scratch, and nothing worth keeping.
    pub tmp: PathBuf,
}

impl Paths {
    /// Where `run` writes: every tool against the set.
    pub fn search(&self) -> PathBuf {
        self.run.join("search")
    }

    /// Where `reject` writes: every tool against the originals of its decoys.
    pub fn reject(&self) -> PathBuf {
        self.run.join("reject")
    }

    /// The label named, with every path resolved against the file that named
    /// it.
    pub fn open(label: &str) -> anyhow::Result<Paths> {
        let file = util::paths::File::open(env!("CARGO_MANIFEST_DIR"))?;
        let p: Paths = file.get(label)?;

        Ok(Paths {
            set: file.at(p.set),
            run: file.at(p.run),
            analysis: file.at(p.analysis),
            figures: file.at(p.figures),
            tmp: file.at(p.tmp),
        })
    }

    /// What to print when no label was named.
    pub fn listing(usage: &str) -> anyhow::Result<String> {
        Ok(util::paths::File::open(env!("CARGO_MANIFEST_DIR"))?.listing(usage))
    }
}

const USAGE: &str = "pid <run|reject|parse|plot> --in <label>";

impl Command {
    fn label(&self) -> Option<&str> {
        match self {
            Command::Run(a) => a.label.as_deref(),
            Command::Reject(a) => a.label.as_deref(),
            Command::Plot(a) => a.label.as_deref(),
            Command::Parse(parse::Cmd::Recall(a)) => a.label.as_deref(),
            Command::Parse(parse::Cmd::Cells(a)) => a.label.as_deref(),
            Command::Parse(parse::Cmd::Score(a)) => a.label.as_deref(),
        }
    }
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    let label = match cli.command.label() {
        Some(label) => label.to_string(),
        None => {
            println!("{}", Paths::listing(USAGE)?);
            return Ok(());
        }
    };
    let paths = Paths::open(&label)?;

    match cli.command {
        Command::Run(args) => run::main(args, &paths),
        Command::Reject(args) => reject::main(args, &paths),
        Command::Parse(cmd) => parse::main(cmd, &paths),
        Command::Plot(args) => plot::main(args, &paths),
    }
}

/// One profmark unit's files, which is what every subcommand opens.
pub struct Inputs {
    pub query_hmm: PathBuf,
    pub query_fa: PathBuf,
    pub query_sto: PathBuf,
    /// One aligned fasta per family. psiblast takes an alignment at a time and
    /// will not read stockholm.
    pub afa: PathBuf,
    pub target_fa: PathBuf,
    /// Which pair is true, and at what identity.
    pub truth: PathBuf,
    /// The decoys unreversed, under the names `target_fa` gives them.
    pub originals: PathBuf,
}

impl Inputs {
    pub fn open(set_dir: &Path) -> anyhow::Result<Inputs> {
        let set = Set::load_as(set_dir, &shape::PROFMARK)?;
        let unit = set
            .units()
            .next()
            .with_context(|| format!("{} names no units", set_dir.display()))?;

        let query_hmm = unit.query_hmm()?;

        Ok(Inputs {
            query_fa: unit.query_fa()?,
            query_sto: unit.query_sto()?,
            // beside the profiles rather than in the manifest: the recipe
            // writes it and no search asks for it
            afa: query_hmm
                .parent()
                .unwrap_or_else(|| Path::new("."))
                .join("afa"),
            target_fa: unit.target()?,
            truth: set_dir.join(unit.need("truth")?),
            originals: set_dir.join(unit.need("originals")?),
            query_hmm,
        })
    }
}
