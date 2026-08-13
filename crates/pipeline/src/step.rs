use crate::cmd::Cmd;
use crate::label;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum OnError {
    #[default]
    Stop,
    Continue,
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
        }
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

    pub(crate) fn halts(&self) -> bool {
        self.on_error == OnError::Stop && self.cmds.iter().any(|c| c.status().failed())
    }
}

impl From<Cmd> for Step {
    fn from(cmd: Cmd) -> Step {
        Step::serial([cmd])
    }
}
