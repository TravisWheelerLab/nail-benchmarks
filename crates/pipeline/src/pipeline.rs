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
        let cores = Cores::read();

        for (s, step) in steps.iter_mut().enumerate() {
            let s_idx = s + 1;
            step.index = Some(s_idx);
            let step_part = label::filename(s_idx, step.name.as_deref());
            let step_label = step.label();
            let inherited = step.cores;

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

                // resolved once here, so nothing downstream has to know a
                // command's core count can come from the step holding it
                if cmd.cores.is_none() {
                    cmd.cores = inherited;
                }

                // the only ask the machine could never satisfy: anything else is
                // a matter of waiting, since commands hand their cores back
                let want = cmd.cores.unwrap_or(0);
                if want > cores.len() {
                    bail!(
                        "{step_label}.{} wants {want} cores, and the machine has {}",
                        cmd.label(),
                        cores.len()
                    );
                }
            }
        }

        Ok(Pipeline {
            steps,
            sinks,
            stderr_dir: if stderr_wanted { maybe_dir } else { None },
            cores,
        })
    }
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
            let mut held = VecDeque::new();

            for cmd in step.cmds() {
                if held.len() >= step.width() {
                    held.pop_front();
                }
                let lease = self.cores.try_acquire(cmd.cores.unwrap_or(0));
                let cpus = lease.as_ref().map(|l| l.cpus()).unwrap_or_default();
                println!("{} {}", cmd.label(), cmd.line_on(cpus));

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

                // a command the batch was cancelled out from under comes back
                // never having run. leaving it alone lets skip_rest announce it
                // once, as skipped, alongside the ones no worker ever picked up
                if !matches!(step.cmds()[k].status(), Status::NotRun) {
                    // giving up on the run is not a reason to leave a batch of
                    // searches running behind us, so this cancels before it
                    // propagates
                    if let Err(e) = sinks.record(step, &step.cmds()[k]) {
                        batch.cancel(scope);
                        return Err(e);
                    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cpu::Cores;
    use crate::step::OnError;
    use std::sync::{Arc, Mutex};

    #[derive(Default)]
    struct Log {
        /// One entry per `record`, so a command announced twice shows up twice.
        records: Vec<(String, String, Status)>,
        finished: usize,
    }

    impl Log {
        fn of(&self, cmd: &str) -> Vec<&Status> {
            self.records
                .iter()
                .filter(|(_, label, _)| label == cmd)
                .map(|(_, _, status)| status)
                .collect()
        }
    }

    /// A sink that writes down everything it is told, and can be asked to fail
    /// when a particular command arrives.
    #[derive(Clone, Default)]
    struct Recorder {
        log: Arc<Mutex<Log>>,
        broken_by: Option<String>,
    }

    impl Recorder {
        fn broken_by(cmd: &str) -> Recorder {
            Recorder {
                broken_by: Some(cmd.to_string()),
                ..Recorder::default()
            }
        }
    }

    impl Sink for Recorder {
        fn record(&mut self, step: &Step, cmd: &Cmd) -> anyhow::Result<()> {
            self.log.lock().unwrap().records.push((
                step.label(),
                cmd.label(),
                cmd.status().clone(),
            ));
            if self.broken_by.as_deref() == Some(cmd.label().as_str()) {
                bail!("the sink cannot write");
            }
            Ok(())
        }

        fn finish(&mut self) -> anyhow::Result<()> {
            self.log.lock().unwrap().finished += 1;
            Ok(())
        }
    }

    fn sh(name: &str, script: &str) -> Cmd {
        Cmd::new("/bin/sh").name(name).arg("-c", script)
    }

    #[test]
    fn a_command_takes_its_core_count_from_the_step_unless_it_has_its_own() {
        let pipeline = PipelineBuilder::new()
            .step(Step::serial([Cmd::new("/a"), Cmd::new("/b").cores(2)]).cores(4))
            .step(Step::serial([Cmd::new("/c")]))
            .no_stderr()
            .build()
            .unwrap();

        assert_eq!(pipeline.steps[0].cmds()[0].cores, Some(4));
        assert_eq!(pipeline.steps[0].cmds()[1].cores, Some(2));
        assert_eq!(pipeline.steps[1].cmds()[0].cores, None);
    }

    #[test]
    fn asking_for_more_cores_than_the_machine_has_fails_the_build() {
        let too_many = Cores::read().len() + 1;
        let built = PipelineBuilder::new()
            .step(Step::serial([Cmd::new("/a").name("greedy")]).name("big").cores(too_many))
            .no_stderr()
            .build();

        let Err(error) = built else {
            panic!("a pipeline that can never be placed should not build");
        };
        let error = error.to_string();

        assert!(error.contains("greedy"), "{error}");
        assert!(error.contains("big"), "{error}");
        assert!(error.contains(&too_many.to_string()), "{error}");
    }

    #[test]
    fn steps_are_numbered_from_one() {
        let pipeline = PipelineBuilder::new()
            .step(Cmd::new("/a"))
            .step(Cmd::new("/b"))
            .no_stderr()
            .build()
            .unwrap();

        assert_eq!(pipeline.steps[0].label(), "[1]");
        assert_eq!(pipeline.steps[1].label(), "[2]");
    }

    #[test]
    fn an_unrouted_stderr_gets_a_file_of_its_own_and_a_routed_one_is_left_alone() {
        let pipeline = PipelineBuilder::new()
            .step(
                Step::serial([
                    Cmd::new("/a").name("x"),
                    Cmd::new("/b").name("y").stderr(Output::Inherit),
                ])
                .name("s"),
            )
            .stderr_dir("/tmp/pipeline-test-stderr")
            .build()
            .unwrap();

        let Output::OnFailure(path) = &pipeline.steps[0].cmds()[0].stderr else {
            panic!("expected a failure log, got {:?}", pipeline.steps[0].cmds()[0].stderr);
        };
        assert!(path.ends_with("1-s.1-x.stderr"), "{}", path.display());
        assert_eq!(pipeline.steps[0].cmds()[1].stderr, Output::Inherit);
        assert!(pipeline.stderr_dir.is_some());
    }

    #[test]
    fn no_directory_is_claimed_when_nothing_needed_a_failure_log() {
        let pipeline = PipelineBuilder::new()
            .step(Cmd::new("/a").stderr(Output::Inherit))
            .stderr_dir("/tmp/pipeline-test-unused")
            .build()
            .unwrap();

        assert!(pipeline.stderr_dir.is_none());
    }

    #[test]
    fn a_failed_step_ends_the_run_and_says_which_command_did_it() {
        let recorder = Recorder::default();
        let log = Arc::clone(&recorder.log);

        let error = PipelineBuilder::new()
            .step(Step::serial([sh("first", "exit 0")]))
            .step(Step::serial([sh("bad", "exit 1")]))
            .step(Step::serial([sh("never", "exit 0")]))
            .no_stderr()
            .sink(recorder)
            .build()
            .unwrap()
            .run()
            .unwrap_err()
            .to_string();

        assert!(error.contains("bad"), "{error}");

        let log = log.lock().unwrap();
        assert!(matches!(log.of("first")[..], [Status::Finished(_)]));
        assert!(matches!(log.of("never")[..], [Status::Skipped]));
    }

    #[test]
    fn continue_lets_the_rest_of_the_run_happen() {
        let recorder = Recorder::default();
        let log = Arc::clone(&recorder.log);

        PipelineBuilder::new()
            .step(Step::serial([sh("bad", "exit 1"), sh("after", "exit 0")]).on_error(OnError::Continue))
            .step(Step::serial([sh("later", "exit 0")]))
            .no_stderr()
            .sink(recorder)
            .build()
            .unwrap()
            .run()
            .expect("a tolerated failure should not end the run");

        let log = log.lock().unwrap();
        assert!(matches!(log.of("after")[..], [Status::Finished(_)]));
        assert!(matches!(log.of("later")[..], [Status::Finished(_)]));
    }

    #[test]
    fn skip_leaves_the_rest_of_its_own_step_and_nothing_else() {
        let recorder = Recorder::default();
        let log = Arc::clone(&recorder.log);

        PipelineBuilder::new()
            .step(
                Step::serial([sh("bad", "exit 1"), sh("sibling", "exit 0")])
                    .on_error(OnError::Skip),
            )
            .step(Step::serial([sh("later", "exit 0")]))
            .no_stderr()
            .sink(recorder)
            .build()
            .unwrap()
            .run()
            .expect("skip should not fail the run");

        let log = log.lock().unwrap();
        assert!(matches!(log.of("sibling")[..], [Status::Skipped]));
        assert!(matches!(log.of("later")[..], [Status::Finished(_)]));
    }

    #[test]
    fn a_sink_that_cannot_write_still_gets_told_the_run_is_over() {
        let recorder = Recorder::broken_by("bad");
        let log = Arc::clone(&recorder.log);

        let error = PipelineBuilder::new()
            .step(Step::serial([sh("bad", "exit 0")]))
            .step(Step::serial([sh("later", "exit 0")]))
            .no_stderr()
            .sink(recorder)
            .build()
            .unwrap()
            .run()
            .unwrap_err()
            .to_string();

        assert!(error.contains("cannot write"), "{error}");
        assert_eq!(
            log.lock().unwrap().finished,
            1,
            "finish must run even when the walk gave up"
        );
    }

    #[test]
    fn every_command_is_announced_exactly_once_and_never_as_not_run() {
        // the contract Sink promises: by the time you see a command it has
        // finished, failed to start, or been skipped
        let recorder = Recorder::default();
        let log = Arc::clone(&recorder.log);

        let names = [
            "serial-ok",
            "batched-a",
            "batched-b",
            "batched-c",
            "stopper",
            "sibling",
            "unreached",
        ];

        PipelineBuilder::new()
            .step(Step::serial([sh("serial-ok", "exit 0")]))
            .step(Step::batched(
                2,
                [
                    sh("batched-a", "exit 0"),
                    sh("batched-b", "exit 0"),
                    sh("batched-c", "exit 0"),
                ],
            ))
            .step(
                Step::batched(1, [sh("stopper", "exit 1"), sh("sibling", "sleep 5")])
                    .on_error(OnError::Skip),
            )
            .step(Step::serial([sh("unreached", "exit 0")]))
            .no_stderr()
            .sink(recorder)
            .build()
            .unwrap()
            .run()
            .expect("skip should not fail the run");

        let log = log.lock().unwrap();
        for name in names {
            let seen = log.of(name);
            assert_eq!(seen.len(), 1, "{name} was announced {} times", seen.len());
            assert!(
                !matches!(seen[0], Status::NotRun),
                "{name} was announced as NotRun"
            );
        }
        assert_eq!(log.records.len(), names.len(), "something extra was announced");
    }

    #[test]
    fn a_command_the_batch_gave_up_on_before_starting_it_is_announced_once() {
        // the awkward case for the contract: a worker takes a command, parks
        // waiting for cores that the hog is holding, and is woken by the cancel.
        // it comes back never having run, and both the batch and skip_rest have
        // an opinion about saying so
        let whole_machine = Cores::read().len();
        let recorder = Recorder::default();
        let log = Arc::clone(&recorder.log);

        PipelineBuilder::new()
            .step(
                Step::batched(
                    3,
                    [
                        // slow enough that the other two are placed first
                        sh("boom", "sleep 0.3; exit 1"),
                        sh("hog", "sleep 5").cores(whole_machine),
                        sh("waiter", "sleep 5").cores(1),
                    ],
                )
                .on_error(OnError::Skip),
            )
            .no_stderr()
            .sink(recorder)
            .build()
            .unwrap()
            .run()
            .expect("skip should not fail the run");

        let log = log.lock().unwrap();
        let seen = log.of("waiter");
        assert_eq!(seen.len(), 1, "announced {} times: {seen:?}", seen.len());
        assert!(!matches!(seen[0], Status::NotRun), "announced as NotRun");
    }

    #[test]
    fn a_batch_reports_every_command_it_was_given() {
        let recorder = Recorder::default();
        let log = Arc::clone(&recorder.log);

        PipelineBuilder::new()
            .step(Step::batched(
                3,
                (0..6).map(|i| sh(&format!("job-{i}"), "exit 0")),
            ))
            .no_stderr()
            .sink(recorder)
            .build()
            .unwrap()
            .run()
            .unwrap();

        let log = log.lock().unwrap();
        assert_eq!(log.records.len(), 6);
        for i in 0..6 {
            assert!(matches!(log.of(&format!("job-{i}"))[..], [Status::Finished(_)]));
        }
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

