//! The list of steps, and running them.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;
use std::time::Instant;

use crate::cmd::Cmd;
use crate::execute::{Status, execute};
use crate::sink::Sink;
use crate::step::{Step, Strategy};

/// A list of steps and somewhere to send what happens.
#[derive(Default)]
pub struct Pipeline {
    steps: Vec<Step>,
    sinks: Sinks,
}

impl Pipeline {
    pub fn new() -> Self {
        Pipeline::default()
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

    /// Print what [`run`](Self::run) would run, without running it. Shows the
    /// exact argv, so a typo in a generated flag is visible before anything
    /// takes an hour to fail.
    pub fn dry_run(&self) {
        for step in &self.steps {
            println!("# {}", step.name());
            for cmd in step.cmds() {
                println!("{}", cmd.line());
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
        } = self;

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

        sinks.finish()
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
    let plan: Vec<Cmd> = step.cmds().to_vec();
    let next = AtomicUsize::new(0);
    let (tx, rx) = mpsc::channel();

    std::thread::scope(|scope| -> anyhow::Result<()> {
        for _ in 0..jobs {
            let tx = tx.clone();
            let next = &next;
            let plan = &plan;
            scope.spawn(move || {
                loop {
                    let k = next.fetch_add(1, Ordering::Relaxed);
                    if k >= plan.len() {
                        break;
                    }
                    let mut cmd = plan[k].clone();
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
