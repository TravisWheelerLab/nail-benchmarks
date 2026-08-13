//! Where results go.
//!
//! A [`Pipeline`](crate::Pipeline) knows nothing about output. It runs
//! commands and announces what happened to whatever sinks were registered on
//! it. Printing progress and writing the summary table are both just sinks.
//!
//! Two things hold for every sink, and they are what make one easy to write:
//! **each command is announced exactly once**, and **never with
//! [`Status::NotRun`]** — by the time you see it, it has finished, failed to
//! start, or been skipped.
//!
//! [`Status::NotRun`]: crate::Status

use crate::cmd::Cmd;
use crate::step::Step;

/// Something that wants to hear what the pipeline did.
///
/// Every method does nothing by default, so an implementation only writes the
/// ones it cares about. Returning `Err` from any of them stops the run: a sink
/// that cannot write is worth abandoning a long benchmark over.
pub trait Sink {
    /// Everything the pipeline is about to run, before any of it has.
    ///
    /// This is where a sink that needs to know the shape of the whole run — the
    /// full set of label keys, say — works it out.
    fn start(&mut self, steps: &[Step]) -> anyhow::Result<()> {
        let _ = steps;
        Ok(())
    }

    /// One command reached its final state. `step` is the one holding it, for
    /// the name and whether its commands ran together.
    fn record(&mut self, step: &Step, cmd: &Cmd) -> anyhow::Result<()> {
        let _ = (step, cmd);
        Ok(())
    }

    /// Every command in `step` has been announced, and its wall clock is final.
    fn step_done(&mut self, step: &Step) -> anyhow::Result<()> {
        let _ = step;
        Ok(())
    }

    /// The pipeline is over, however it ended.
    fn finish(&mut self) -> anyhow::Result<()> {
        Ok(())
    }
}
