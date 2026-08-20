//! What a run is doing, a line at a time.
//!
//! This is the pipeline's output rather than a note about it, so it goes to
//! stdout. Nothing in this crate writes to stderr.

use std::io::Write;

use crate::execute::Status;
use crate::fmt::{bytes, dash};
use crate::item::Item;
use crate::sink::Sink;
use crate::step::Step;

#[derive(Debug, Default)]
pub struct Progress {
    titled: bool,
}

impl Progress {
    pub fn new() -> Progress {
        Progress::default()
    }

    fn title(&mut self, step: &Step) {
        if !self.titled {
            self.titled = true;
            say(&step.label());
        }
    }
}

/// One line out.
///
/// Flushed rather than left to the buffer, because stdout holds its output back
/// in blocks when it is not a terminal, and a benchmark that runs for an hour
/// is worth being able to watch through a pipe.
fn say(line: &str) {
    let mut out = std::io::stdout().lock();
    let _ = writeln!(out, "{line}");
    let _ = out.flush();
}

/// One line about how something went, whichever kind of thing it was.
fn outcome(name: &str, status: &Status) -> String {
    match (status, status.timing()) {
        (Status::NotRun, _) => format!("  {name:<32} not run"),
        (Status::Skipped, _) => format!("  {name:<32} skipped"),
        (Status::Failed(why), _) => format!("  {name:<32} {why}"),
        (status, Some(t)) => format!(
            "  {:<32} {:>9.2}s {:>9}  {}",
            name,
            t.wall_s,
            t.max_rss_kb.map(bytes).unwrap_or_else(dash),
            match status {
                Status::TimedOut(_) => "timed out".to_string(),
                _ if t.ok() => "ok".to_string(),
                _ => format!("exit {}", t.exit),
            }
        ),
        (_, None) => format!("  {name:<32} -"),
    }
}

impl Sink for Progress {
    fn item_done(&mut self, step: &Step, _at: usize, item: Item<'_>) -> anyhow::Result<()> {
        self.title(step);
        say(&outcome(&item.label(), item.status()));

        // only if it is still there — a failure that said nothing has had its
        // file cleaned up already
        if item.status().failed()
            && let Some(path) = item.stderr_path()
            && path.exists()
        {
            say(&format!("    {}", path.display()));
        }
        Ok(())
    }

    fn step_done(&mut self, _step: &Step) -> anyhow::Result<()> {
        self.titled = false;
        Ok(())
    }
}
