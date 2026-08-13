//! A sink that says on stderr what is happening as it happens.

use crate::cmd::Cmd;
use crate::execute::Status;
use crate::sink::Sink;
use crate::step::Step;
use crate::table::format_bytes;

/// Prints a line per command as it lands, under a header per step.
#[derive(Debug, Default)]
pub struct Progress {
    steps: usize,
    index: usize,
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
    fn start(&mut self, steps: &[Step]) -> anyhow::Result<()> {
        self.steps = steps.len();
        Ok(())
    }

    fn record(&mut self, step: &Step, cmd: &Cmd) -> anyhow::Result<()> {
        if !self.titled {
            self.index += 1;
            self.titled = true;
            eprintln!("[{}/{}] {}", self.index, self.steps, step.name());
        }

        let name = cmd.name();
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
            }
        );
        Ok(())
    }

    fn step_done(&mut self, _step: &Step) -> anyhow::Result<()> {
        self.titled = false;
        Ok(())
    }
}
