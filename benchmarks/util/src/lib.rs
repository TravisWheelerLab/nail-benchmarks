//! What the benchmarks in this repo agree on.
//!
//! Each benchmark asks its own question and keeps its own analyses. What is
//! the same for all of them is the shape of the record each stage leaves
//! behind, and where that record goes.
//!
//! A build describes what it produced in [`set`], one row per search unit. A
//! run is recorded by a `michi::Table` sink as `manifest.tbl`, which [`ledger`]
//! distills into the part an analysis needs. What the analysis works out gets
//! written through [`tbl`], which is the format all three ride on.
//!
//! Nothing here decides where any of that goes. [`paths`] reads the file each
//! tool keeps beside its own source, naming what it reads and what it writes,
//! so a location is something a tool is told rather than something a library
//! knows.
//!
//! [`set`] and [`ledger`] are the same bargain at two seams: a fixed spine plus
//! an open map, written by the producer in the pass that produces the artifact.
//! That is what lets a search read a set without knowing which recipe built it,
//! and an analysis read a run without knowing which pipeline ran it.
//!
//! Where the programs being benchmarked and the sequence data they run on were
//! put is [`tools`].
//!
//! Two smaller things are shared for want of a second home: the pair that cut
//! a query set up, [`split`] by weight, for a batch of jobs, and [`cut`] by
//! name, for a subset or a file per record. [`time`] is
//! for the other way a run can be recorded: timed by a shell rather than by a
//! pipeline, on a machine this workspace never sees.
//!
//! Nothing here decides what counts as a *true* hit. That is the part every
//! benchmark answers differently, and it stays with the benchmark.

pub mod clean;
pub mod cut;
pub mod ledger;
pub mod manifest;
pub mod paths;
pub mod set;
pub mod split;
pub mod tbl;
pub mod time;
pub mod tools;

/// A profile of `leng` nodes over a two-symbol alphabet, which is the smallest
/// thing libsail's parser accepts.
//
// pub because the tests that need a profile libsail will read are in two
// crates: `split`'s here, and mgy's around the Pfam cuts
pub fn profile(name: &str, leng: usize) -> String {
    let mut out = format!(
        "HMMER3/f\nNAME  {name}\nLENG  {leng}\nALPH  amino\n\
         HMM          A        C\n\
         \x20           m->m     m->i     m->d     i->m     i->i     d->m     d->d\n\
         \x20      1.00000  1.00000\n\
         \x20      0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000\n"
    );

    for node in 1..=leng {
        out.push_str(&format!("{node:>7}  1.00000  1.00000\n"));
        out.push_str("       1.00000  1.00000\n");
        out.push_str("       0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000\n");
    }

    out.push_str("//\n");
    out
}
