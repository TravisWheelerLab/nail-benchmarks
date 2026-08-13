//! The list of steps, and running them.

use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;
use std::time::Instant;

use anyhow::Context;

use crate::cmd::{Cmd, Output};
use crate::execute::{Status, execute, stamp};
use crate::label;
use crate::sink::Sink;
use crate::step::{Step, Strategy};

/// Where failure logs go unless told otherwise.
const STDERR_DIR: &str = "stderr";

/// A pipeline under construction.
pub struct PipelineBuilder {
    steps: Vec<Step>,
    sinks: Sinks,
    stderr_dir: Option<PathBuf>,
}

impl Default for PipelineBuilder {
    fn default() -> PipelineBuilder {
        PipelineBuilder {
            steps: Vec::new(),
            sinks: Sinks::default(),
            stderr_dir: Some(PathBuf::from(STDERR_DIR)),
        }
    }
}

impl PipelineBuilder {
    pub fn new() -> Self {
        PipelineBuilder::default()
    }

    /// Add a step. A bare [`Cmd`](crate::Cmd) counts as one, so `.step(cmd)`
    /// works.
    pub fn step(mut self, step: impl Into<Step>) -> Self {
        self.steps.push(step.into());
        self
    }

    /// Add somewhere for results to go. A pipeline with no sinks runs silently
    /// and writes nothing.
    pub fn sink(mut self, sink: impl Sink + 'static) -> Self {
        self.sinks.0.push(Box::new(sink));
        self
    }

    /// Keep failure logs somewhere other than `stderr/`. Each run gets its own
    /// directory under this one, named after the time it was built.
    pub fn stderr_dir(mut self, dir: impl Into<PathBuf>) -> Self {
        self.stderr_dir = Some(dir.into());
        self
    }

    /// Do not keep stderr from failures at all. Commands that route their own
    /// stderr are unaffected either way.
    pub fn no_stderr(mut self) -> Self {
        self.stderr_dir = None;
        self
    }

    /// Settle everything that can be settled before anything runs: number the
    /// steps and their commands, and decide where a failure's stderr would go.
    ///
    /// Nothing here touches the disk, so building and then only looking — see
    /// [`Pipeline::dry_run`] — leaves nothing behind.
    pub fn build(self) -> Pipeline {
        let PipelineBuilder {
            mut steps,
            sinks,
            stderr_dir,
        } = self;

        let dir = stderr_dir.map(|dir| dir.join(stamp()));
        let mut wanted = false;

        for (s, step) in steps.iter_mut().enumerate() {
            let s = s + 1;
            step.index = Some(s);
            let step_part = label::filename(s, step.name.as_deref());

            for (c, cmd) in step.cmds_mut().iter_mut().enumerate() {
                let c = c + 1;
                cmd.index = Some((s, c));

                // a command that has not said where its stderr goes gets a file
                // of its own, which it only keeps if it fails
                if let Some(dir) = &dir
                    && cmd.stderr == Output::Null
                {
                    let cmd_part = label::filename(c, cmd.name.as_deref());
                    cmd.stderr =
                        Output::OnFailure(dir.join(format!("{step_part}.{cmd_part}.stderr")));
                    wanted = true;
                }
            }
        }

        Pipeline {
            steps,
            sinks,
            stderr_dir: if wanted { dir } else { None },
        }
    }
}

/// A pipeline with everything decided, ready to run.
pub struct Pipeline {
    steps: Vec<Step>,
    sinks: Sinks,
    /// This run's directory for failure logs, if any command is going to use it.
    stderr_dir: Option<PathBuf>,
}

impl Pipeline {
    /// Print what [`run`](Self::run) would run, without running it. Shows the
    /// exact argv, so a typo in a generated flag is visible before anything
    /// takes an hour to fail.
    pub fn dry_run(&self) {
        for step in &self.steps {
            println!("# {}", step.label());
            for cmd in step.cmds() {
                println!("{} {}", cmd.label(), cmd.line());
            }
        }
    }

    /// Run everything, announcing each command to every sink as it lands.
    ///
    /// A command failing is not an error: it gets a [`Status`] and the step's
    /// [`OnError`](crate::OnError) decides whether the run goes on. When it does
    /// not, every command left — in this step and every step after it — is
    /// announced as [`Status::Skipped`], so a sink never has to work out what it
    /// missed.
    ///
    /// `Err` means the pipeline could not do its job, which in practice means a
    /// sink refused to write.
    pub fn run(self) -> anyhow::Result<()> {
        let Pipeline {
            mut steps,
            mut sinks,
            stderr_dir,
        } = self;

        // a redirect target has to exist before anything spawns, so this cannot
        // wait for the first failure
        if let Some(dir) = &stderr_dir {
            std::fs::create_dir_all(dir)
                .with_context(|| format!("failed to create {}", dir.display()))?;
        }

        sinks.start(&steps)?;

        // run steps until one of them says that is enough
        let mut left = steps.iter_mut();
        for step in left.by_ref() {
            match step.strategy {
                Strategy::Serial => serial(step, &mut sinks)?,
                Strategy::Batch { jobs } => batch(step, jobs, &mut sinks)?,
            }
            skip_rest(step, &mut sinks)?;
            sinks.step_done(step)?;

            if step.halts() {
                break;
            }
        }

        // the ones that never got a turn still have to be accounted for
        for step in left {
            skip_rest(step, &mut sinks)?;
            sinks.step_done(step)?;
        }

        sinks.finish()?;

        // nothing failed loudly, so there is nothing in there; both of these do
        // nothing if the directory has anything in it
        if let Some(dir) = &stderr_dir {
            std::fs::remove_dir(dir).ok();
            if let Some(parent) = dir.parent() {
                std::fs::remove_dir(parent).ok();
            }
        }
        Ok(())
    }
}

/// Every sink registered on a pipeline, so the rest of this file can talk to
/// them as if there were one.
#[derive(Default)]
struct Sinks(Vec<Box<dyn Sink>>);

impl Sinks {
    fn start(&mut self, steps: &[Step]) -> anyhow::Result<()> {
        self.0.iter_mut().try_for_each(|s| s.start(steps))
    }

    fn record(&mut self, step: &Step, cmd: &Cmd) -> anyhow::Result<()> {
        self.0.iter_mut().try_for_each(|s| s.record(step, cmd))
    }

    fn step_done(&mut self, step: &Step) -> anyhow::Result<()> {
        self.0.iter_mut().try_for_each(|s| s.step_done(step))
    }

    fn finish(&mut self) -> anyhow::Result<()> {
        self.0.iter_mut().try_for_each(|s| s.finish())
    }
}

/// One command at a time, stopping early if the step says to.
fn serial(step: &mut Step, sinks: &mut Sinks) -> anyhow::Result<()> {
    let start = Instant::now();

    for j in 0..step.cmds().len() {
        execute(&mut step.cmds_mut()[j]);
        step.elapsed_s = Some(start.elapsed().as_secs_f64());

        sinks.record(step, &step.cmds()[j])?;
        if step.halts() {
            break;
        }
    }
    Ok(())
}

/// `jobs` at a time. Workers take the next command whenever they are free, so
/// one slow command does not idle the rest.
fn batch(step: &mut Step, jobs: usize, sinks: &mut Sinks) -> anyhow::Result<()> {
    let start = Instant::now();

    // the workers run against a copy, which leaves the step free for the main
    // thread to write results into; cloning a Cmd costs nothing next to
    // spawning a process
    let copies: Vec<Cmd> = step.cmds().to_vec();
    let next = AtomicUsize::new(0);
    let (tx, rx) = mpsc::channel();

    std::thread::scope(|scope| -> anyhow::Result<()> {
        for _ in 0..jobs {
            let tx = tx.clone();
            let next = &next;
            let copies = &copies;
            scope.spawn(move || {
                loop {
                    let k = next.fetch_add(1, Ordering::Relaxed);
                    if k >= copies.len() {
                        break;
                    }
                    let mut cmd = copies[k].clone();
                    execute(&mut cmd);
                    // a closed channel means the main thread gave up on us
                    if tx.send((k, cmd.status)).is_err() {
                        break;
                    }
                }
            });
        }
        drop(tx);

        // results are applied here rather than in the workers, so sinks stay
        // single-threaded and the table keeps up during a long batch
        for (k, status) in rx {
            step.cmds_mut()[k].status = status;
            step.elapsed_s = Some(start.elapsed().as_secs_f64());
            sinks.record(step, &step.cmds()[k])?;
        }
        Ok(())
    })
}

/// Marks anything that never ran as skipped and tells the sinks about it.
///
/// Two cases end up here: the tail of a step that stopped partway, and every
/// command of a step the pipeline never got to.
fn skip_rest(step: &mut Step, sinks: &mut Sinks) -> anyhow::Result<()> {
    for j in 0..step.cmds().len() {
        if matches!(step.cmds()[j].status(), Status::NotRun) {
            step.cmds_mut()[j].status = Status::Skipped;
            sinks.record(step, &step.cmds()[j])?;
        }
    }
    Ok(())
}
