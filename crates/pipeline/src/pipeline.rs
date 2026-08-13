//! The list of steps, and running them.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;
use std::time::Instant;

use crate::cmd::Cmd;
use crate::execute::{execute, Status};
use crate::sink::Sink;
use crate::step::{OnError, Step};

/// A list of steps and somewhere to send what happens.
#[derive(Default)]
pub struct Pipeline {
    steps: Vec<Step>,
    sinks: Vec<Box<dyn Sink>>,
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
        self.sinks.push(Box::new(sink));
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
    /// [`OnError`] decides whether the run goes on. When it does not, every
    /// command left — in this step and every step after it — is announced as
    /// [`Status::Skipped`], so a sink never has to work out what it missed.
    ///
    /// `Err` means the pipeline could not do its job, which in practice means a
    /// sink refused to write.
    pub fn run(self) -> anyhow::Result<()> {
        let Pipeline {
            mut steps,
            mut sinks,
        } = self;

        for sink in &mut sinks {
            sink.start(&steps)?;
        }

        let mut stopped = false;
        for step in &mut steps {
            if !stopped {
                let start = Instant::now();
                match step.jobs() {
                    Some(jobs) => batch(step, jobs, start, &mut sinks)?,
                    None => serial(step, start, &mut sinks)?,
                }
                stopped = step.failed() && step.on_error == OnError::Stop;
            }

            skip_rest(step, &mut sinks)?;
            for sink in &mut sinks {
                sink.step_done(step)?;
            }
        }

        for sink in &mut sinks {
            sink.finish()?;
        }
        Ok(())
    }
}

/// One command at a time, stopping early if the step says to.
fn serial(step: &mut Step, start: Instant, sinks: &mut [Box<dyn Sink>]) -> anyhow::Result<()> {
    for j in 0..step.cmds().len() {
        execute(&mut step.cmds_mut()[j]);
        step.elapsed_s = Some(start.elapsed().as_secs_f64());

        let cmd = &step.cmds()[j];
        for sink in sinks.iter_mut() {
            sink.record(step, cmd)?;
        }

        if step.cmds()[j].status().failed() && step.on_error == OnError::Stop {
            break;
        }
    }
    Ok(())
}

/// `jobs` at a time. Workers take the next command whenever they are free, so
/// one slow command does not idle the rest.
fn batch(
    step: &mut Step,
    jobs: usize,
    start: Instant,
    sinks: &mut [Box<dyn Sink>],
) -> anyhow::Result<()> {
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

            let cmd = &step.cmds()[k];
            for sink in sinks.iter_mut() {
                sink.record(step, cmd)?;
            }
        }
        Ok(())
    })
}

/// Announce whatever never got its turn, so every command is heard from once.
fn skip_rest(step: &mut Step, sinks: &mut [Box<dyn Sink>]) -> anyhow::Result<()> {
    let mut skipped = Vec::new();
    for j in 0..step.cmds().len() {
        if matches!(step.cmds()[j].status(), Status::NotRun) {
            step.cmds_mut()[j].status = Status::Skipped;
            skipped.push(j);
        }
    }

    for j in skipped {
        let cmd = &step.cmds()[j];
        for sink in sinks.iter_mut() {
            sink.record(step, cmd)?;
        }
    }
    Ok(())
}
