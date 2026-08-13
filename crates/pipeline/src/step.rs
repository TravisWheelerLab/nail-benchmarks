use crate::cmd::Cmd;
use crate::label;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum OnError {
    #[default]
    Stop,
    Continue,
}

/// How a step's commands are run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strategy {
    /// One after another, in the order given.
    Serial,
    /// `jobs` at a time, in no particular order.
    Batch { jobs: usize },
}

#[derive(Clone, Debug)]
pub struct Step {
    pub(crate) name: Option<String>,
    /// Where this step sits in the pipeline, counting from one. Handed out by
    /// [`PipelineBuilder::build`](crate::PipelineBuilder::build).
    pub(crate) index: Option<usize>,
    cmds: Vec<Cmd>,
    pub(crate) strategy: Strategy,
    pub(crate) on_error: OnError,
    /// Wall clock around the whole step, set as it runs. `None` until it
    /// starts. This is measured, not added up from the commands: it catches
    /// spawn overhead and the gaps between them, and for a batch it is the real
    /// elapsed time however many waves the jobs took.
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

    /// What to call this step in output. Optional: without one it goes by its
    /// number alone.
    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    pub fn on_error(mut self, on_error: OnError) -> Self {
        self.on_error = on_error;
        self
    }

    /// How this step is named in output: `[2](search)`, or `[2]` if it was never
    /// given a name.
    pub fn label(&self) -> String {
        let index = self.index.map(|index| index.to_string());
        label::label(index.as_deref(), self.name.as_deref())
    }

    pub fn cmds(&self) -> &[Cmd] {
        &self.cmds
    }

    /// How long the step took on the wall clock, or `None` if it never ran.
    pub fn wall_s(&self) -> Option<f64> {
        self.elapsed_s
    }

    /// Whether the commands ran one after another or together, and how many at
    /// a time if together.
    pub fn strategy(&self) -> Strategy {
        self.strategy
    }

    pub(crate) fn cmds_mut(&mut self) -> &mut [Cmd] {
        &mut self.cmds
    }

    /// Whether something went wrong here and this step's policy says that ends
    /// the run. The same question a serial step asks after every command and the
    /// pipeline asks after every step.
    pub(crate) fn halts(&self) -> bool {
        self.on_error == OnError::Stop && self.cmds.iter().any(|c| c.status().failed())
    }
}

/// A bare command is a step of one. The step stays unnamed, so the two do not
/// print the same name twice on the row they collapse into.
impl From<Cmd> for Step {
    fn from(cmd: Cmd) -> Step {
        Step::serial([cmd])
    }
}
