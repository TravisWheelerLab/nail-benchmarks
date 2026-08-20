use crate::cmd::Cmd;
use crate::label;

/// How far a failed command reaches.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum OnError {
    /// Run the rest of this step's commands anyway.
    Continue,
    /// Skip the rest of this step. The pipeline carries on to the next one.
    Skip,
    /// Skip the rest of this step, skip every step after it, and fail the run.
    #[default]
    Abort,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strategy {
    Serial,
    Batched { jobs: usize },
}

#[derive(Clone, Debug)]
pub struct Step {
    pub(crate) name: Option<String>,
    pub(crate) index: Option<usize>,
    cmds: Vec<Cmd>,
    pub(crate) strategy: Strategy,
    pub(crate) on_error: OnError,
    pub(crate) elapsed_s: Option<f64>,
    /// How many cores each of this step's commands asks for, unless it asked
    /// for itself.
    pub(crate) cores: Option<usize>,
}

impl Step {
    pub fn serial(cmds: impl IntoIterator<Item = Cmd>) -> Self {
        Step::new(cmds, Strategy::Serial)
    }

    pub fn batched(jobs: usize, cmds: impl IntoIterator<Item = Cmd>) -> Self {
        Step::new(cmds, Strategy::Batched { jobs: jobs.max(1) })
    }

    fn new(cmds: impl IntoIterator<Item = Cmd>, strategy: Strategy) -> Self {
        Step {
            name: None,
            index: None,
            cmds: cmds.into_iter().collect(),
            strategy,
            on_error: OnError::default(),
            elapsed_s: None,
            cores: None,
        }
    }

    /// Pin each of this step's commands to `cores` physical cores.
    ///
    /// Per command, not per step: a batch of four with `cores(2)` asks for eight
    /// cores. If the machine cannot spare that many at once, the commands that
    /// cannot be placed wait for the ones that can.
    pub fn cores(mut self, cores: usize) -> Self {
        self.cores = Some(cores);
        self
    }

    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    pub fn on_error(mut self, on_error: OnError) -> Self {
        self.on_error = on_error;
        self
    }

    /// What to call this step. Empty until the pipeline is built, which is when
    /// a step learns its number.
    pub fn label(&self) -> String {
        match self.index {
            Some(index) => label::label(index, self.name.as_deref()),
            None => self.name.clone().unwrap_or_default(),
        }
    }

    pub fn cmds(&self) -> &[Cmd] {
        &self.cmds
    }

    pub fn wall_s(&self) -> Option<f64> {
        self.elapsed_s
    }

    pub fn strategy(&self) -> Strategy {
        self.strategy
    }

    /// How many of this step's commands can be running at once.
    pub fn width(&self) -> usize {
        match self.strategy {
            Strategy::Serial => 1,
            Strategy::Batched { jobs } => jobs,
        }
    }

    pub(crate) fn cmds_mut(&mut self) -> &mut [Cmd] {
        &mut self.cmds
    }

    /// Whether the rest of this step's commands are worth running.
    pub(crate) fn skips(&self) -> bool {
        self.on_error != OnError::Continue && self.cmds.iter().any(|c| c.status().failed())
    }

    /// The command that ends the run, if this step holds one.
    pub(crate) fn aborts(&self) -> Option<&Cmd> {
        (self.on_error == OnError::Abort)
            .then(|| self.cmds.iter().find(|c| c.status().failed()))
            .flatten()
    }
}

impl From<Cmd> for Step {
    fn from(cmd: Cmd) -> Step {
        Step::serial([cmd])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::execute::{Status, Timing};

    fn timing(exit: i32) -> Timing {
        Timing {
            wall_s: 1.0,
            user_s: 1.0,
            sys_s: 0.0,
            max_rss_kb: 1024,
            exit,
        }
    }

    /// Two commands, the second having gone however `outcome` says.
    fn step(on_error: OnError, outcome: Status) -> Step {
        let mut step = Step::serial([Cmd::new("/a"), Cmd::new("/b")]).on_error(on_error);
        step.cmds[0].status = Status::Finished(timing(0));
        step.cmds[1].status = outcome;
        step
    }

    #[test]
    fn nothing_failed_so_nothing_reaches_anywhere() {
        for on_error in [OnError::Continue, OnError::Skip, OnError::Abort] {
            let step = step(on_error, Status::Finished(timing(0)));
            assert!(!step.skips(), "{on_error:?} skipped a clean step");
            assert!(step.aborts().is_none(), "{on_error:?} aborted a clean step");
        }
    }

    #[test]
    fn continue_lets_a_failure_pass() {
        let step = step(OnError::Continue, Status::Finished(timing(1)));
        assert!(!step.skips());
        assert!(step.aborts().is_none());
    }

    #[test]
    fn skip_stops_the_step_and_stops_there() {
        let step = step(OnError::Skip, Status::Finished(timing(1)));
        assert!(step.skips());
        assert!(step.aborts().is_none());
    }

    #[test]
    fn abort_stops_the_step_and_names_what_did_it() {
        let step = step(OnError::Abort, Status::Finished(timing(1)));
        assert!(step.skips());
        assert_eq!(step.aborts().map(|cmd| cmd.label()), Some("b".to_string()));
    }

    #[test]
    fn abort_is_the_default() {
        assert_eq!(OnError::default(), OnError::Abort);
    }

    #[test]
    fn a_command_that_never_ran_is_not_a_failure() {
        for outcome in [Status::NotRun, Status::Skipped] {
            let step = step(OnError::Abort, outcome);
            assert!(!step.skips());
            assert!(step.aborts().is_none());
        }
    }

    #[test]
    fn every_way_of_failing_counts() {
        for outcome in [
            Status::Failed("could not spawn".into()),
            Status::TimedOut(timing(143)),
            Status::Finished(timing(1)),
        ] {
            let step = step(OnError::Abort, outcome);
            assert!(step.skips());
            assert!(step.aborts().is_some());
        }
    }

    #[test]
    fn a_batch_always_has_a_worker() {
        // zero jobs would mean no workers at all, which hangs rather than
        // finishing empty
        assert_eq!(Step::batched(0, [Cmd::new("/a")]).width(), 1);
        assert_eq!(Step::batched(4, [Cmd::new("/a")]).width(), 4);
        assert_eq!(Step::serial([Cmd::new("/a")]).width(), 1);
    }
}
