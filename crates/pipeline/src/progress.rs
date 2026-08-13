//! A sink that says on stderr what is happening as it happens.

use crate::cmd::Cmd;
use crate::execute::Status;
use crate::sink::Sink;
use crate::step::Step;
use crate::table::format_bytes;

/// Prints a line per command as it lands, under a header per step.
#[derive(Debug, Default)]
pub struct Progress {
    /// Whether the current step's header has been printed. The header waits for
    /// the step's first command so it needs no hook of its own.
    titled: bool,
}

impl Progress {
    pub fn new() -> Progress {
        Progress::default()
    }
}

impl Sink for Progress {
    fn record(&mut self, step: &Step, cmd: &Cmd) -> anyhow::Result<()> {
        if !self.titled {
            self.titled = true;
            eprintln!("{}", step.label());
        }

        let name = cmd.label();
        eprintln!(
            "{}",
            match cmd.status() {
                Status::NotRun => format!("  {name:<32} not run"),
                Status::Skipped => format!("  {name:<32} skipped"),
                Status::Failed(why) => format!("  {name:<32} {why}"),
                Status::Finished(t) => format!(
                    "  {:<32} {:>9.2}s {:>9}  {}",
                    name,
                    t.wall_s,
                    format_bytes(t.max_rss_kb),
                    if t.ok() {
                        "ok".to_string()
                    } else {
                        format!("exit {}", t.exit)
                    }
                ),
                Status::TimedOut(t) => format!(
                    "  {:<32} {:>9.2}s {:>9}  timed out",
                    name,
                    t.wall_s,
                    format_bytes(t.max_rss_kb)
                ),
            }
        );

        // only if it is still there — a failure that said nothing has had its
        // file cleaned up already
        if cmd.status().failed()
            && let Some(path) = cmd.stderr_path()
            && path.exists()
        {
            eprintln!("    {}", path.display());
        }
        Ok(())
    }

    fn step_done(&mut self, _step: &Step) -> anyhow::Result<()> {
        self.titled = false;
        Ok(())
    }
}
