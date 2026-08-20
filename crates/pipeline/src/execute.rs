//! Running commands and finding out what they cost.

use std::fs::{File, OpenOptions};
use std::io;
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Condvar, Mutex};
use std::thread::Scope;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use anyhow::Context;

use crate::cmd::{Cmd, Output};
use crate::cpu::{Cores, Lease};

/// How long anything gets to leave on a SIGTERM before it gets a SIGKILL.
const SIGTERM_GRACE: Duration = Duration::from_secs(5);

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Timing {
    pub wall_s: f64,
    pub user_s: f64,
    pub sys_s: f64,
    pub max_rss_kb: i64,
    pub exit: i32,
}

impl Timing {
    pub fn ok(&self) -> bool {
        self.exit == 0
    }
}

#[derive(Clone, Debug)]
pub enum Status {
    NotRun,
    Skipped,
    Failed(String),
    Finished(Timing),
    TimedOut(Timing),
}

impl Status {
    pub fn failed(&self) -> bool {
        match self {
            Status::NotRun | Status::Skipped => false,
            Status::Failed(_) | Status::TimedOut(_) => true,
            Status::Finished(t) => !t.ok(),
        }
    }

    /// What this cost, for the two states that got far enough to have an
    /// answer. A command killed on its deadline still burned everything it says
    /// it burned, so it counts.
    pub fn timing(&self) -> Option<&Timing> {
        match self {
            Status::Finished(t) | Status::TimedOut(t) => Some(t),
            _ => None,
        }
    }
}

impl Cores {
    /// Run one command. Nothing can cut it short but its own timeout.
    pub(crate) fn execute(&self, cmd: &mut Cmd) {
        // held until the command is finished with, however it finishes; the cpus
        // go back on the way out of scope
        let lease = self.acquire(cmd.cores.unwrap_or(0), &|| false);
        cmd.cpus = pinned(&lease);

        let status = match cmd.spawn() {
            Err(e) => Status::Failed(format!("{e:#}")),
            Ok((pid, start)) => outcome(wait(pid, start, cmd.timeout, || {})),
        };
        cmd.report(status);
    }
}

/// One batch step: whether it has been cancelled, and the pids it has running.
///
/// Both are made fresh for the step, so cancelling reaches that step's commands
/// and nothing else, and leaves nothing behind for the next step.
pub(crate) struct Batch<'a> {
    cores: &'a Cores,
    cancelled: AtomicBool,
    /// A cancel has to reach every process at once, and this is the only place
    /// their pids are collected.
    running: Mutex<Vec<libc::pid_t>>,
    emptied: Condvar,
}

impl<'a> Batch<'a> {
    pub(crate) fn new(cores: &'a Cores) -> Batch<'a> {
        Batch {
            cores,
            cancelled: AtomicBool::new(false),
            running: Mutex::new(Vec::new()),
            emptied: Condvar::new(),
        }
    }

    pub(crate) fn cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed)
    }

    /// Run one command as part of the batch, which may be cancelled under it.
    pub(crate) fn execute(&self, cmd: &mut Cmd) {
        let lease = self.cores.acquire(cmd.cores.unwrap_or(0), &|| self.cancelled());

        // an empty lease means one of two things: the command asked for no
        // pinning, or the wait was cut short by a cancel. only the second is a
        // reason not to run it — and leaving the status alone reports it as
        // skipped, like the commands no worker reached
        if self.cancelled() {
            return;
        }
        cmd.cpus = pinned(&lease);

        let status = match cmd.spawn() {
            Err(e) => Status::Failed(format!("{e:#}")),
            Ok((pid, start)) => {
                self.add(pid);
                let waited = wait(pid, start, cmd.timeout, || self.remove(pid));
                // wait removes it on the way through, but not if it failed first
                self.remove(pid);
                outcome(waited)
            }
        };
        cmd.report(status);
    }

    /// Nothing new starts, and everything running is terminated. Calling this
    /// more than once does nothing extra, which a batch relies on: it cancels
    /// once per result after the first failure.
    pub(crate) fn cancel<'s>(&'s self, scope: &'s Scope<'s, '_>) {
        if self.cancelled.swap(true, Ordering::Relaxed) {
            return;
        }

        // the flag is set before the lock, so a command registering right now
        // either appears in this list or sees the flag itself
        let running = self.running.lock().unwrap();
        for pid in running.iter() {
            signal(*pid, libc::SIGTERM);
        }
        drop(running);

        // a worker asleep waiting for cores would never see the flag on its own,
        // and the cores it is waiting for may never come back
        self.cores.wake();

        scope.spawn(|| self.kill_remaining());
    }

    /// Whatever ignores the SIGTERM gets a SIGKILL once the grace is up. Only
    /// pids still listed, and a listed pid has not been reaped, so this cannot
    /// reach a process that stopped being ours.
    fn kill_remaining(&self) {
        let running = self.running.lock().unwrap();
        let (running, grace) = self
            .emptied
            .wait_timeout_while(running, SIGTERM_GRACE, |running| !running.is_empty())
            .unwrap();

        if grace.timed_out() {
            for pid in running.iter() {
                signal(*pid, libc::SIGKILL);
            }
        }
    }

    fn add(&self, pid: libc::pid_t) {
        let mut running = self.running.lock().unwrap();
        running.push(pid);
        // one that started while the batch was being cancelled still has to go
        if self.cancelled() {
            signal(pid, libc::SIGTERM);
        }
    }

    fn remove(&self, pid: libc::pid_t) {
        let mut running = self.running.lock().unwrap();
        running.retain(|p| *p != pid);
        if running.is_empty() {
            self.emptied.notify_all();
        }
    }
}

fn pinned(lease: &Option<Lease>) -> Vec<usize> {
    match lease {
        Some(lease) => lease.cpus().to_vec(),
        None => Vec::new(),
    }
}

fn outcome(waited: anyhow::Result<(Timing, bool)>) -> Status {
    match waited {
        Ok((timing, true)) => Status::TimedOut(timing),
        Ok((timing, false)) => Status::Finished(timing),
        Err(e) => Status::Failed(format!("{e:#}")),
    }
}

fn signal(pid: libc::pid_t, sig: libc::c_int) {
    unsafe { libc::kill(pid, sig) };
}

pub(crate) fn stamp() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    let now = secs as libc::time_t;
    let mut tm: libc::tm = unsafe { std::mem::zeroed() };
    let mut buf = [0u8; 32];

    let written = unsafe {
        if libc::localtime_r(&now, &mut tm).is_null() {
            return secs.to_string();
        }
        libc::strftime(
            buf.as_mut_ptr().cast(),
            buf.len(),
            c"%Y%m%d-%H%M%S".as_ptr(),
            &tm,
        )
    };
    String::from_utf8_lossy(&buf[..written]).into_owned()
}

impl Cmd {
    /// Start the process, giving back its pid and the moment it started.
    fn spawn(&self) -> anyhow::Result<(libc::pid_t, Instant)> {
        let (program, args) = self.argv();
        let mut proc = Command::new(&program);
        proc.args(args);
        proc.envs(&self.env);
        if let Some(dir) = &self.dir {
            proc.current_dir(dir);
        }
        proc.stdout(self.stdout.stdio()?);
        proc.stderr(self.stderr.stdio()?);

        let start = Instant::now();
        let child = proc
            .spawn()
            .with_context(|| format!("failed to spawn {}", self.program.display()))?;

        Ok((child.id() as libc::pid_t, start))
    }

    /// Record how the command went, and drop a stderr file that was only being
    /// kept in case of failure and has nothing in it.
    fn report(&mut self, status: Status) {
        if let Output::OnFailure(path) = &self.stderr {
            let keep = status.failed() && std::fs::metadata(path).is_ok_and(|m| m.len() > 0);
            if !keep {
                std::fs::remove_file(path).ok();
            }
        }
        self.status = status;
    }
}

impl Output {
    fn stdio(&self) -> anyhow::Result<Stdio> {
        let make_dir = |path: &Path| -> anyhow::Result<()> {
            if let Some(dir) = path.parent() {
                std::fs::create_dir_all(dir)?;
            }
            Ok(())
        };

        Ok(match self {
            Output::Null => Stdio::null(),
            Output::Inherit => Stdio::inherit(),
            Output::File(path) | Output::OnFailure(path) => {
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
}

/// A process with a deadline on it.
///
/// Only made for a command that asked for a timeout — one that did not is waited
/// on directly and never touches a lock. `finished` is set while the process is
/// over but not yet reaped, which is the window in which the waiting thread can
/// be told to leave the pid alone.
struct Deadline {
    pid: libc::pid_t,
    finished: Mutex<bool>,
    changed: Condvar,
}

impl Deadline {
    fn new(pid: libc::pid_t) -> Deadline {
        Deadline {
            pid,
            finished: Mutex::new(false),
            changed: Condvar::new(),
        }
    }

    /// Wait `after` for the process, then terminate it. True if it came to that.
    fn kill_after(&self, after: Duration) -> bool {
        if self.wait_for_finish(after) {
            return false;
        }
        signal(self.pid, libc::SIGTERM);
        if !self.wait_for_finish(SIGTERM_GRACE) {
            signal(self.pid, libc::SIGKILL);
        }
        true
    }

    /// Wait up to `limit` for the process to be over. True if it is.
    fn wait_for_finish(&self, limit: Duration) -> bool {
        let finished = self.finished.lock().unwrap();
        let (finished, _) = self
            .changed
            .wait_timeout_while(finished, limit, |finished| !*finished)
            .unwrap();
        *finished
    }

    fn set_finished(&self) {
        *self.finished.lock().unwrap() = true;
        self.changed.notify_all();
    }
}

/// Wait for the process, terminating it if it runs longer than `limit`. The flag
/// says whether it came to that.
///
/// `exited` is handed straight to [`reap`]. A command with no limit needs no
/// second thread and no lock, which is the common case.
fn wait(
    pid: libc::pid_t,
    start: Instant,
    limit: Option<Duration>,
    exited: impl FnOnce(),
) -> anyhow::Result<(Timing, bool)> {
    let Some(limit) = limit else {
        return reap(pid, start, exited).map(|timing| (timing, false));
    };

    let deadline = Deadline::new(pid);
    std::thread::scope(|scope| {
        let waiting = scope.spawn(|| deadline.kill_after(limit));
        let timing = reap(pid, start, || {
            deadline.set_finished();
            exited();
        });
        // reap sets it on the way through, but not if it failed before it got
        // there, and the other thread would sit on the limit for nothing
        deadline.set_finished();
        let killed = waiting.join().unwrap_or(false);
        timing.map(|timing| (timing, killed))
    })
}

/// Block until the process is over, then collect what it cost.
///
/// Two waits rather than one. `waitid` says it is over but leaves the pid held,
/// and only `wait4` hands it back — `exited` runs in the gap between them, which
/// is the one moment anything else holding this pid can be told to stop
/// signalling it.
fn reap(pid: libc::pid_t, start: Instant, exited: impl FnOnce()) -> anyhow::Result<Timing> {
    let mut info: libc::siginfo_t = unsafe { std::mem::zeroed() };
    let rc = unsafe {
        libc::waitid(
            libc::P_PID,
            pid as libc::id_t,
            &mut info,
            libc::WEXITED | libc::WNOWAIT,
        )
    };
    if rc < 0 {
        return Err(io::Error::last_os_error()).context("waitid failed");
    }
    // taken here rather than after the reap, so it is when the process ended
    let wall_s = start.elapsed().as_secs_f64();

    exited();

    let mut status: libc::c_int = 0;
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    // std never reaps a child on its own, so nothing in it races with this
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
        wall_s,
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
