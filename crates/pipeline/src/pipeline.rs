//! The list of steps, and running them.

use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, mpsc};
use std::time::Instant;

use anyhow::Context;

use crate::cmd::{Cmd, Output};
use crate::execute::{Live, Status, stamp};
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

    pub fn step(mut self, step: impl Into<Step>) -> Self {
        self.steps.push(step.into());
        self
    }

    pub fn sink(mut self, sink: impl Sink + 'static) -> Self {
        self.sinks.0.push(Box::new(sink));
        self
    }

    pub fn stderr_dir(mut self, dir: impl Into<PathBuf>) -> Self {
        self.stderr_dir = Some(dir.into());
        self
    }

    pub fn no_stderr(mut self) -> Self {
        self.stderr_dir = None;
        self
    }

    pub fn build(self) -> Pipeline {
        let PipelineBuilder {
            mut steps,
            sinks,
            stderr_dir,
        } = self;

        let maybe_dir = stderr_dir.map(|dir| dir.join(stamp()));
        let mut stderr_wanted = false;

        for (s, step) in steps.iter_mut().enumerate() {
            let s_idx = s + 1;
            step.index = Some(s_idx);
            let step_part = label::filename(s_idx, step.name.as_deref());

            for (c, cmd) in step.cmds_mut().iter_mut().enumerate() {
                let c_idx = c + 1;
                cmd.index = Some((s_idx, c_idx));

                // if we have a stderr dir, we redirect every un-routed stderr
                // to its own file just in case the command ends up failing
                if let Some(dir) = &maybe_dir
                    && cmd.stderr == Output::Null
                {
                    let cmd_part = label::filename(c_idx, cmd.name.as_deref());
                    cmd.stderr =
                        Output::OnFailure(dir.join(format!("{step_part}.{cmd_part}.stderr")));
                    stderr_wanted = true;
                }
            }
        }

        Pipeline {
            steps,
            sinks,
            stderr_dir: if stderr_wanted { maybe_dir } else { None },
        }
    }
}

pub struct Pipeline {
    steps: Vec<Step>,
    sinks: Sinks,
    stderr_dir: Option<PathBuf>,
}

impl Pipeline {
    pub fn dry_run(&self) {
        for step in &self.steps {
            println!("# {}", step.label());
            for cmd in step.cmds() {
                println!("{} {}", cmd.label(), cmd.line());
            }
        }
    }

    pub fn run(self) -> anyhow::Result<()> {
        let Pipeline {
            mut steps,
            mut sinks,
            stderr_dir,
        } = self;

        sinks.start(&steps)?;

        if let Some(dir) = &stderr_dir {
            std::fs::create_dir_all(dir)
                .with_context(|| format!("failed to create {}", dir.display()))?;
        }

        let live = Arc::new(Live::default());
        let _watchdog = Live::watch(&live);

        // run steps until one of them halts
        let mut left = steps.iter_mut();
        for step in left.by_ref() {
            match step.strategy {
                Strategy::Serial => serial(step, &live, &mut sinks)?,
                Strategy::Batch { jobs } => batch(step, jobs, &live, &mut sinks)?,
            }
            skip_rest(step, &mut sinks)?;
            sinks.step_done(step)?;

            if step.halts() {
                break;
            }
        }

        // any steps after an early halt get to report here
        for step in left {
            skip_rest(step, &mut sinks)?;
            sinks.step_done(step)?;
        }

        sinks.finish()?;

        // cleanup
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
fn serial(step: &mut Step, live: &Live, sinks: &mut Sinks) -> anyhow::Result<()> {
    let start = Instant::now();

    for j in 0..step.cmds().len() {
        live.execute(&mut step.cmds_mut()[j]);
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
///
/// A step that has failed and is set to stop does not wait for the rest of the
/// batch to finish: nothing new is taken, and what is already running is killed.
/// Those come back as `exit 143`, which is a real failure and reads as one, so a
/// step that stopped shows the one command that broke it and the ones it took
/// down with it.
fn batch(step: &mut Step, jobs: usize, live: &Live, sinks: &mut Sinks) -> anyhow::Result<()> {
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
                    if live.stopping() {
                        break;
                    }
                    let k = next.fetch_add(1, Ordering::Relaxed);
                    if k >= copies.len() {
                        break;
                    }
                    let mut cmd = copies[k].clone();
                    live.execute(&mut cmd);
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

            if step.halts() {
                live.stop();
            }
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
