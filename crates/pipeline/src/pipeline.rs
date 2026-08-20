//! The list of steps, and running them.

use std::collections::VecDeque;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;
use std::time::Instant;

use anyhow::{Context, anyhow, bail};

use crate::cmd::{Cmd, Output};
use crate::cpu::Cores;
use crate::execute::{Batch, Status, stamp};
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

    pub fn build(self) -> anyhow::Result<Pipeline> {
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

        // resolved once here, so nothing downstream has to know a command's
        // core count can come from the step holding it
        for step in steps.iter_mut() {
            let inherited = step.cores;
            for cmd in step.cmds_mut() {
                if cmd.cores.is_none() {
                    cmd.cores = inherited;
                }
            }
        }

        let cores = Cores::read();
        fits(&steps, &cores)?;

        Ok(Pipeline {
            steps,
            sinks,
            stderr_dir: if stderr_wanted { maybe_dir } else { None },
            cores,
        })
    }
}

/// Whether every command asks for something the machine could ever give it.
///
/// Only the one check: a command wanting more cores than exist would wait for
/// them forever. Anything else is a matter of waiting, since commands hand
/// their cores back when they finish.
fn fits(steps: &[Step], cores: &Cores) -> anyhow::Result<()> {
    for step in steps {
        for cmd in step.cmds() {
            let want = cmd.cores.unwrap_or(0);
            if want > cores.len() {
                bail!(
                    "{}.{} wants {want} cores, and the machine has {}",
                    step.label(),
                    cmd.label(),
                    cores.len()
                );
            }
        }
    }
    Ok(())
}

pub struct Pipeline {
    steps: Vec<Step>,
    sinks: Sinks,
    stderr_dir: Option<PathBuf>,
    cores: Cores,
}

impl Pipeline {
    pub fn dry_run(&self) {
        for step in &self.steps {
            println!("# {}", step.label());

            // cores really are taken and given back here, so the pinning shown
            // is one a run could produce. only as many are held at once as the
            // step could run at once; which command gets which is a guess,
            // since that depends on the order workers pick them up
            let width = match step.strategy {
                Strategy::Serial => 1,
                Strategy::Batched { jobs } => jobs,
            };
            let mut held = VecDeque::new();

            for cmd in step.cmds() {
                if held.len() >= width {
                    held.pop_front();
                }
                let lease = self.cores.try_acquire(cmd.cores.unwrap_or(0));

                let mut copy = cmd.clone();
                copy.cpus = lease.as_ref().map(|l| l.cpus().to_vec()).unwrap_or_default();
                println!("{} {}", copy.label(), copy.line());

                held.extend(lease);
            }
        }
    }

    pub fn run(mut self) -> anyhow::Result<()> {
        // the steps come out so that the rest of the pipeline can be borrowed
        // while they are worked through. run consumes self, so nobody sees the
        // field it leaves behind
        let mut steps = std::mem::take(&mut self.steps);

        self.sinks.start(&steps)?;

        if let Some(dir) = &self.stderr_dir {
            std::fs::create_dir_all(dir)
                .with_context(|| format!("failed to create {}", dir.display()))?;
        }

        let outcome = self.run_steps(&mut steps);

        // both of these happen whichever way the run went, so the table is
        // finished and the stderr log is where it says it is by the time the
        // caller hears about anything
        let finished = self.sinks.finish();
        if let Some(dir) = &self.stderr_dir {
            std::fs::remove_dir(dir).ok();
            if let Some(parent) = dir.parent() {
                std::fs::remove_dir(parent).ok();
            }
        }

        // a sink that could not write comes first: if the table never made it to
        // disk, which command failed is not on record anyway
        let failure = outcome?;
        finished?;

        // last, so everything above has already happened
        match failure {
            Some(failure) => Err(anyhow!(failure)),
            None => Ok(()),
        }
    }

    /// Every step in turn, stopping at the first one that ends the run. Gives
    /// back the command that ended it, if one did.
    fn run_steps(&mut self, steps: &mut [Step]) -> anyhow::Result<Option<String>> {
        let mut failure = None;
        let mut remaining_steps = steps.iter_mut();

        for step in remaining_steps.by_ref() {
            match step.strategy {
                Strategy::Serial => self.serial(step)?,
                Strategy::Batched { jobs } => self.batch(step, jobs)?,
            }
            self.skip_rest(step)?;
            self.sinks.step_done(step)?;

            if let Some(cmd) = step.aborts() {
                failure = Some(format!("{} failed: {}", cmd.label(), cmd.line()));
                break;
            }
        }

        // any steps after an early abort get to report here
        for step in remaining_steps {
            self.skip_rest(step)?;
            self.sinks.step_done(step)?;
        }

        Ok(failure)
    }

    /// One command at a time, stopping early if the step says to.
    fn serial(&mut self, step: &mut Step) -> anyhow::Result<()> {
        let start = Instant::now();

        for j in 0..step.cmds().len() {
            self.cores.execute(&mut step.cmds_mut()[j]);
            step.elapsed_s = Some(start.elapsed().as_secs_f64());

            self.sinks.record(step, &step.cmds()[j])?;
            if step.skips() {
                break;
            }
        }
        Ok(())
    }

    /// `jobs` at a time. Workers take the next command whenever they are free, so
    /// one slow command does not idle the rest.
    ///
    /// A step that has failed and is set to stop does not wait for the rest of
    /// the batch to finish: nothing new is taken, and what is already running is
    /// killed. Those come back as `exit 143`, which is a real failure and reads
    /// as one, so a step that stopped shows the one command that broke it and the
    /// ones it took down with it.
    fn batch(&mut self, step: &mut Step, jobs: usize) -> anyhow::Result<()> {
        let start = Instant::now();

        // the workers run against a copy, which leaves the step free for the main
        // thread to write results into; cloning a Cmd costs nothing next to
        // spawning a process
        let copies: Vec<Cmd> = step.cmds().to_vec();
        let next = AtomicUsize::new(0);
        let (tx, rx) = mpsc::channel();
        // made here and dropped here, so a cancel reaches this step's commands
        // and nothing else
        let batch = Batch::new(&self.cores);
        let sinks = &mut self.sinks;

        std::thread::scope(|scope| -> anyhow::Result<()> {
            for _ in 0..jobs {
                let tx = tx.clone();
                let next = &next;
                let copies = &copies;
                let batch = &batch;
                scope.spawn(move || {
                    loop {
                        if batch.cancelled() {
                            break;
                        }
                        let k = next.fetch_add(1, Ordering::Relaxed);
                        if k >= copies.len() {
                            break;
                        }
                        let mut cmd = copies[k].clone();
                        batch.execute(&mut cmd);
                        // a closed channel means the main thread gave up on us
                        if tx.send((k, cmd)).is_err() {
                            break;
                        }
                    }
                });
            }
            drop(tx);

            // results are applied here rather than in the workers, so sinks stay
            // single-threaded and the table keeps up during a long batch
            // the whole command comes back, not just its status, so the cpus the
            // worker actually pinned it to are what the table reports
            for (k, cmd) in rx {
                step.cmds_mut()[k] = cmd;
                step.elapsed_s = Some(start.elapsed().as_secs_f64());

                // giving up on the run is not a reason to leave a batch of
                // searches running behind us, so this cancels before it propagates
                if let Err(e) = sinks.record(step, &step.cmds()[k]) {
                    batch.cancel(scope);
                    return Err(e);
                }

                if step.skips() {
                    batch.cancel(scope);
                }
            }
            Ok(())
        })
    }

    /// Marks anything that never ran as skipped and tells the sinks about it.
    ///
    /// Two cases end up here: the tail of a step that stopped partway, and every
    /// command of a step the pipeline never got to.
    fn skip_rest(&mut self, step: &mut Step) -> anyhow::Result<()> {
        for j in 0..step.cmds().len() {
            if matches!(step.cmds()[j].status(), Status::NotRun) {
                step.cmds_mut()[j].status = Status::Skipped;
                self.sinks.record(step, &step.cmds()[j])?;
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

