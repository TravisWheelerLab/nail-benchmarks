use std::path::PathBuf;
use std::time::Instant;

use crate::execute::{execute_batch, execute_verbose};
use crate::step::{OnError, Step};
use crate::table::Report;

#[derive(Clone, Debug, Default)]
pub struct Pipeline {
    steps: Vec<Step>,
    table: Option<PathBuf>,
}

impl Pipeline {
    pub fn new() -> Self {
        Pipeline::default()
    }

    pub fn step(mut self, step: impl Into<Step>) -> Self {
        self.steps.push(step.into());
        self
    }

    pub fn table(mut self, path: impl Into<PathBuf>) -> Self {
        self.table = Some(path.into());
        self
    }

    pub fn dry_run(&self) {
        for step in &self.steps {
            println!("# {}", step.name());
            for cmd in step.cmds() {
                println!("{}", cmd.line());
            }
        }
    }

    pub fn run(self) -> anyhow::Result<Report> {
        let mut report = Report { steps: self.steps };
        let total = report.steps.len();

        if let Some(path) = &self.table {
            report.write(path)?;
        }

        'pipeline: for i in 0..total {
            eprintln!("[{}/{}] {}", i + 1, total, report.steps[i].name());
            let start = Instant::now();

            match report.steps[i].jobs() {
                Some(jobs) => {
                    execute_batch(report.steps[i].cmds_mut(), jobs);
                    report.steps[i].elapsed_s = Some(start.elapsed().as_secs_f64());

                    if let Some(path) = &self.table {
                        report.write(path)?;
                    }
                    if report.steps[i].failed() && report.steps[i].on_error == OnError::Stop {
                        eprintln!("  stopping: step {}", report.steps[i].name());
                        break 'pipeline;
                    }
                }
                None => {
                    for j in 0..report.steps[i].cmds().len() {
                        let cmd = &mut report.steps[i].cmds_mut()[j];
                        execute_verbose(cmd);
                        let failed = cmd.status().failed();
                        // updated per command, so a step still in progress
                        // shows how long it has been going
                        report.steps[i].elapsed_s = Some(start.elapsed().as_secs_f64());

                        if let Some(path) = &self.table {
                            report.write(path)?;
                        }

                        if failed && report.steps[i].on_error == OnError::Stop {
                            eprintln!("  stopping: step {}", report.steps[i].name());
                            break 'pipeline;
                        }
                    }
                }
            }
        }

        Ok(report)
    }
}
