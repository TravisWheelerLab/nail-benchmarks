use crate::closure::Closure;
use crate::cmd::Cmd;
use crate::execute::Status;
use crate::fmt::{bytes, dash};
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
            eprintln!("{}", step.label());
        }
    }
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
    fn record(&mut self, step: &Step, cmd: &Cmd) -> anyhow::Result<()> {
        self.title(step);
        eprintln!("{}", outcome(&cmd.label(), cmd.status()));

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

    fn record_closure(&mut self, step: &Step, closure: &Closure) -> anyhow::Result<()> {
        self.title(step);
        eprintln!("{}", outcome(closure.label(), closure.status()));
        Ok(())
    }

    fn step_done(&mut self, _step: &Step) -> anyhow::Result<()> {
        self.titled = false;
        Ok(())
    }
}
