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
    Batch { jobs: usize },
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

    pub fn batch(jobs: usize, cmds: impl IntoIterator<Item = Cmd>) -> Self {
        Step::new(cmds, Strategy::Batch { jobs: jobs.max(1) })
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
    /// Per command, not per step: a batch of four with `cores(2)` has eight
    /// cores in flight, on eight distinct physical cores.
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

    pub fn label(&self) -> String {
        let index = self.index.map(|index| index.to_string());
        label::label(index.as_deref(), self.name.as_deref())
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
