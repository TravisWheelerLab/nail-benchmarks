//! Running commands and finding out what they cost.

use std::fs::{File, OpenOptions};
use std::io;
use std::process::{Command, Stdio};
use std::time::Instant;

use anyhow::Context;

use crate::cmd::{Cmd, Output};

/// What one command cost, taken from `wait4`.
///
/// The numbers come out of the kernel rather than being formatted to text by
/// something like `/usr/bin/time` and parsed back.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Timing {
    pub wall_s: f64,
    pub user_s: f64,
    pub sys_s: f64,
    pub max_rss_kb: i64,
    /// The process's exit code, or 128 + signal if it was killed.
    pub exit: i32,
}

impl Timing {
    pub fn ok(&self) -> bool {
        self.exit == 0
    }
}

/// How far a command got.
#[derive(Clone, Debug)]
pub enum Status {
    /// Not yet attempted. The starting state, and the only one a [`Sink`] never
    /// sees: by the time a command is announced it has reached one of the
    /// three below.
    ///
    /// [`Sink`]: crate::Sink
    NotRun,
    /// Never will be attempted, because the pipeline stopped before reaching it.
    Skipped,
    /// Tried, but no process came of it — the program could not be spawned, or
    /// a redirect could not be opened. Carries the reason.
    Failed(String),
    /// A process ran and was measured. It may still have exited nonzero; that
    /// is in [`Timing::exit`].
    Finished(Timing),
}

impl Status {
    /// Whether this counts against the run. A command that never got its turn
    /// did not fail.
    pub fn failed(&self) -> bool {
        match self {
            Status::NotRun | Status::Skipped => false,
            Status::Failed(_) => true,
            Status::Finished(t) => !t.ok(),
        }
    }
}

/// Run one command to completion and hang its [`Status`] on it.
///
/// Nothing here is an error. A nonzero exit comes back in [`Timing::exit`], and
/// a command that could never start comes back as [`Status::Failed`], so a
/// failure to launch is a row in the table rather than a gap.
///
/// Running an already-run command just overwrites what was there.
pub fn execute(cmd: &mut Cmd) {
    let status = match spawn_and_wait(cmd) {
        Ok(timing) => Status::Finished(timing),
        Err(e) => Status::Failed(format!("{e:#}")),
    };
    cmd.status = status;
}

fn spawn_and_wait(cmd: &Cmd) -> anyhow::Result<Timing> {
    let mut proc = Command::new(&cmd.argv[0]);
    proc.args(&cmd.argv[1..]);
    proc.stdout(stdio(&cmd.stdout)?);
    proc.stderr(stdio(&cmd.stderr)?);

    let start = Instant::now();
    let child = proc
        .spawn()
        .with_context(|| format!("failed to spawn {}", cmd.argv[0]))?;

    wait_timed(child.id() as libc::pid_t, start)
}

/// Open a redirect target, making its directory if it is not there yet.
fn stdio(out: &Output) -> anyhow::Result<Stdio> {
    let make_dir = |path: &std::path::Path| -> anyhow::Result<()> {
        if let Some(dir) = path.parent() {
            std::fs::create_dir_all(dir)?;
        }
        Ok(())
    };

    Ok(match out {
        Output::Null => Stdio::null(),
        Output::Inherit => Stdio::inherit(),
        Output::File(path) => {
            make_dir(path)?;
            let file = File::create(path)
                .with_context(|| format!("failed to create {}", path.display()))?;
            Stdio::from(file)
        }
        Output::Append(path) => {
            make_dir(path)?;
            let file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)
                .with_context(|| format!("failed to open {} for append", path.display()))?;
            Stdio::from(file)
        }
    })
}

/// Reap this specific child with `wait4` so its rusage comes back on its own.
/// `getrusage(RUSAGE_CHILDREN)` would fold in every sibling, which stops being
/// right the moment anything runs more than one process.
fn wait_timed(pid: libc::pid_t, start: Instant) -> anyhow::Result<Timing> {
    let mut status: libc::c_int = 0;
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };

    // std::process::Child does not reap on drop, so taking the child here does
    // not race with anything in std.
    let rc = unsafe { libc::wait4(pid, &mut status, 0, &mut usage) };
    if rc < 0 {
        return Err(io::Error::last_os_error()).context("wait4 failed");
    }

    let exit = if libc::WIFEXITED(status) {
        libc::WEXITSTATUS(status)
    } else if libc::WIFSIGNALED(status) {
        128 + libc::WTERMSIG(status)
    } else {
        -1
    };

    Ok(Timing {
        wall_s: start.elapsed().as_secs_f64(),
        user_s: secs(usage.ru_utime),
        sys_s: secs(usage.ru_stime),
        // ru_maxrss is already in kilobytes on Linux
        max_rss_kb: usage.ru_maxrss,
        exit,
    })
}

fn secs(tv: libc::timeval) -> f64 {
    tv.tv_sec as f64 + tv.tv_usec as f64 / 1_000_000.0
}
