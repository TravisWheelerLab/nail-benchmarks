use crate::cmd::Cmd;
use crate::execute::Status;
use crate::fmt::bytes;
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
            match (cmd.status(), cmd.status().timing()) {
                (Status::NotRun, _) => format!("  {name:<32} not run"),
                (Status::Skipped, _) => format!("  {name:<32} skipped"),
                (Status::Failed(why), _) => format!("  {name:<32} {why}"),
                (status, Some(t)) => format!(
                    "  {:<32} {:>9.2}s {:>9}  {}",
                    name,
                    t.wall_s,
                    bytes(t.max_rss_kb),
                    match status {
                        Status::TimedOut(_) => "timed out".to_string(),
                        _ if t.ok() => "ok".to_string(),
                        _ => format!("exit {}", t.exit),
                    }
                ),
                (_, None) => format!("  {name:<32} -"),
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
