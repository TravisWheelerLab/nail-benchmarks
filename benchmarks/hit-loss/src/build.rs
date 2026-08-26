//! Draws the one query set and the one target set a run measures against.
//!
//! Copied from cloud-search's build: a single fixed pair of files, not a
//! ladder. This benchmark isn't moving anything about the search's size or
//! shape, so there is nothing to sweep here either.

use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, bail};
use clap::Parser;

use bioio::aggregate::AggregateFasta;
use pail::{Closure, Cmd, PipelineBuilder, Progress, Step};
use tools::{mgnify, pfam_hmm};

pub const DEFAULT_NAME: &str = "benchmark";

// the index that comes with the collection. building one from scratch is a pass
// over every byte of MGnify, so this is not a file to go without
const INDEX_NAME: &str = "mgnify.afi";

#[derive(Parser, Debug)]
pub struct Args {
    /// The name of the benchmark directory, resolved under benchmarks/hit-loss/.
    #[arg(default_value = DEFAULT_NAME)]
    pub name: String,

    /// How many Pfam families to take, off the front of the file. Leave it off
    /// for all of them.
    #[arg(long, value_name = "N")]
    pub families: Option<usize>,

    /// How many MGnify sequences to draw.
    #[arg(long, default_value_t = 100_000, value_name = "N")]
    pub seqs: usize,

    /// Random seed
    #[arg(long, default_value_t = 67779, value_name = "N")]
    pub seed: u64,
}

pub fn main(args: Args) -> anyhow::Result<()> {
    let src_hmm = pfam_hmm()?;
    let src_dir = mgnify()?;

    let bench = crate::dir().join(&args.name);

    if bench.exists() {
        bail!("benchmark: {} already exists", args.name)
    }

    let queries = bench.join("queries");
    let targets = bench.join("targets");
    let query_hmm = queries.join("query.hmm");
    let target_fa = targets.join("target.fa");

    PipelineBuilder::new()
        .step(
            Cmd::new("mkdir")
                .name("dirs")
                .flag("-p")
                .path(&queries)
                .path(&targets),
        )
        .step(
            Step::from_closures([
                Closure::new("query", {
                    let families = args.families;

                    move || {
                        match families {
                            // all of Pfam, so there is nothing to pick out, and
                            // a copy beats reading 1.6GB a line at a time
                            None => {
                                std::fs::copy(&src_hmm, &query_hmm).with_context(|| {
                                    format!("failed to copy {}", src_hmm.display())
                                })?;
                            }
                            Some(n) => {
                                bioio::hmm::subset(&src_hmm, n, &query_hmm)?;
                            }
                        }

                        Ok(())
                    }
                }),
                Closure::new("target", {
                    let (n, seed) = (args.seqs, args.seed);

                    move || {
                        // allow_overwrite only bites when the index no longer
                        // matches its sources, and then it rebuilds rather than
                        // refusing. that is a pass over every byte of MGnify,
                        // so a run that suddenly goes quiet for hours is this
                        let seqs = AggregateFasta::builder()
                            .dir(&src_dir)
                            .index(src_dir.join(INDEX_NAME))
                            .allow_overwrite()
                            .build()?;

                        let mut out = BufWriter::new(File::create(&target_fa)?);
                        seqs.sample(n, seed, &mut out)?;
                        out.flush()?;

                        Ok(())
                    }
                }),
            ])
            .name("draw"),
        )
        .stderr_dir(bench.join("stderr"))
        .sink(Progress::new())
        .build()?
        .run()?;

    println!("\nbuilt {}", bench.display());
    Ok(())
}
