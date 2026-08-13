use crate::cmd::Cmd;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum OnError {
    #[default]
    Stop,
    Continue,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Strategy {
    Serial,
    Batch { jobs: usize },
}

#[derive(Clone, Debug)]
pub struct Step {
    name: String,
    cmds: Vec<Cmd>,
    strategy: Strategy,
    pub(crate) on_error: OnError,
    /// Wall clock around the whole step, set as it runs. `None` until it
    /// starts. This is measured, not added up from the commands: it catches
    /// spawn overhead and the gaps between them, and for a batch it is the real
    /// elapsed time however many waves the jobs took.
    pub(crate) elapsed_s: Option<f64>,
}

impl Step {
    pub fn serial(name: impl Into<String>, cmds: impl IntoIterator<Item = Cmd>) -> Self {
        Step {
            name: name.into(),
            on_error: OnError::default(),
            cmds: cmds.into_iter().collect(),
            strategy: Strategy::Serial,
            elapsed_s: None,
        }
    }

    pub fn batch(
        name: impl Into<String>,
        jobs: usize,
        cmds: impl IntoIterator<Item = Cmd>,
    ) -> Self {
        Step {
            name: name.into(),
            on_error: OnError::default(),
            cmds: cmds.into_iter().collect(),
            strategy: Strategy::Batch { jobs: jobs.max(1) },
            elapsed_s: None,
        }
    }

    pub fn on_error(mut self, on_error: OnError) -> Self {
        self.on_error = on_error;
        self
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn cmds(&self) -> &[Cmd] {
        &self.cmds
    }

    /// How long the step took on the wall clock, or `None` if it never ran.
    pub fn wall_s(&self) -> Option<f64> {
        self.elapsed_s
    }

    /// Whether the commands ran together rather than one after another.
    pub(crate) fn batched(&self) -> bool {
        matches!(self.strategy, Strategy::Batch { .. })
    }

    pub(crate) fn cmds_mut(&mut self) -> &mut [Cmd] {
        &mut self.cmds
    }

    pub(crate) fn jobs(&self) -> Option<usize> {
        match self.strategy {
            Strategy::Serial => None,
            Strategy::Batch { jobs } => Some(jobs),
        }
    }

    pub(crate) fn failed(&self) -> bool {
        self.cmds.iter().any(|c| c.status().failed())
    }
}

impl From<Cmd> for Step {
    fn from(cmd: Cmd) -> Step {
        Step::serial(cmd.name().to_string(), [cmd])
    }
}
