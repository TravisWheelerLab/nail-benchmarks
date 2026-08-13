//! Running commands and finding out what they cost.

use std::fs::{File, OpenOptions};
use std::io;
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use anyhow::Context;

use crate::cmd::{Cmd, Output};

/// How often the watchdog looks at the children.
const TICK: Duration = Duration::from_millis(100);
/// How long a command gets to leave on its own before it is made to.
const GRACE: Duration = Duration::from_secs(5);

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
    /// Ran past its [`Cmd::timeout`] and was killed for it. Everything measured
    /// up to the kill is real, so a command that timed out still says what it
    /// was costing when it went.
    TimedOut(Timing),
}

impl Status {
    /// Whether this counts against the run. A command that never got its turn
    /// did not fail.
    pub fn failed(&self) -> bool {
        match self {
            Status::NotRun | Status::Skipped => false,
            Status::Failed(_) | Status::TimedOut(_) => true,
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
    let live = Arc::new(Live::default());
    // nothing to watch for unless there is a deadline: on its own, one command
    // has nobody who might call it off
    let _watchdog = cmd.timeout.is_some().then(|| Live::watch(&live));
    live.execute(cmd);
}

/// The children running right now, so somebody who is not waiting on them can
/// still end them.
///
/// Two things want that: a step that has already failed and should not keep
/// paying for the rest of its batch, and a command that has run past its
/// deadline. Both come down to signalling a pid somebody else is blocked on.
#[derive(Default)]
pub(crate) struct Live {
    running: Mutex<Vec<Entry>>,
    stopping: AtomicBool,
    done: AtomicBool,
}

/// One running child, and how far along the way out it is.
struct Entry {
    pid: libc::pid_t,
    /// When it has overstayed, for a command that set a timeout.
    deadline: Option<Instant>,
    phase: Phase,
    /// Whether what ended it was its own deadline rather than the pipeline
    /// giving up, which is the difference between `time` and `fail`.
    expired: bool,
}

/// How hard we have asked a child to leave.
enum Phase {
    Running,
    /// Sent SIGTERM at this point, and waiting out [`GRACE`].
    Term(Instant),
    /// Sent SIGKILL. There is nothing after this.
    Kill,
}

impl Live {
    /// Like [`execute`], with the child visible to [`stop`](Self::stop) and to
    /// the watchdog while it runs.
    pub(crate) fn execute(&self, cmd: &mut Cmd) {
        let status = match spawn(cmd) {
            Err(e) => Status::Failed(format!("{e:#}")),
            Ok((pid, start)) => {
                self.add(pid, cmd.timeout.map(|after| start + after));
                let waited = wait_timed(pid, start);
                let expired = self.remove(pid);
                match waited {
                    Ok(t) if expired => Status::TimedOut(t),
                    Ok(t) => Status::Finished(t),
                    Err(e) => Status::Failed(format!("{e:#}")),
                }
            }
        };
        if let Output::OnFailure(path) = &cmd.stderr {
            tidy(path, &status);
        }
        cmd.status = status;
    }

    /// Nothing new starts, and anything already running is asked to leave.
    pub(crate) fn stop(&self) {
        // the lock is held across the flag and the signals, so a child spawning
        // right now either shows up in this list or sees the flag in `add`
        let mut running = self.running.lock().unwrap();
        self.stopping.store(true, Ordering::Relaxed);
        let now = Instant::now();
        for entry in running.iter_mut() {
            entry.term(now);
        }
    }

    pub(crate) fn stopping(&self) -> bool {
        self.stopping.load(Ordering::Relaxed)
    }

    /// Start the thread that enforces deadlines and finishes off anything that
    /// ignored a polite SIGTERM. It stops when the returned handle drops.
    pub(crate) fn watch(live: &Arc<Live>) -> Watchdog {
        let live = Arc::clone(live);
        let handle = std::thread::spawn({
            let live = Arc::clone(&live);
            move || {
                while !live.done.load(Ordering::Relaxed) {
                    live.sweep();
                    std::thread::sleep(TICK);
                }
            }
        });
        Watchdog {
            live,
            handle: Some(handle),
        }
    }

    fn add(&self, pid: libc::pid_t, deadline: Option<Instant>) {
        let mut running = self.running.lock().unwrap();
        let mut entry = Entry {
            pid,
            deadline,
            phase: Phase::Running,
            expired: false,
        };
        if self.stopping() {
            entry.term(Instant::now());
        }
        running.push(entry);
    }

    /// Take a finished child off the list, saying whether its deadline is what
    /// ended it.
    fn remove(&self, pid: libc::pid_t) -> bool {
        let mut running = self.running.lock().unwrap();
        match running.iter().position(|e| e.pid == pid) {
            Some(i) => running.swap_remove(i).expired,
            None => false,
        }
    }

    /// One pass over the children: anything past its deadline gets asked to
    /// leave, and anything that was asked long enough ago stops being asked.
    fn sweep(&self) {
        let now = Instant::now();
        let mut running = self.running.lock().unwrap();
        for entry in running.iter_mut() {
            match entry.phase {
                Phase::Running if entry.deadline.is_some_and(|at| now >= at) => {
                    entry.expired = true;
                    entry.term(now);
                }
                Phase::Term(at) if now.duration_since(at) >= GRACE => {
                    entry.phase = Phase::Kill;
                    signal(entry.pid, libc::SIGKILL);
                }
                _ => {}
            }
        }
    }
}

impl Entry {
    fn term(&mut self, now: Instant) {
        if matches!(self.phase, Phase::Running) {
            self.phase = Phase::Term(now);
            signal(self.pid, libc::SIGTERM);
        }
    }
}

/// Keeps the watchdog thread alive. Dropping it stops the thread and waits for
/// it, so it cannot outlive the run that started it.
pub(crate) struct Watchdog {
    live: Arc<Live>,
    handle: Option<JoinHandle<()>>,
}

impl Drop for Watchdog {
    fn drop(&mut self) {
        self.live.done.store(true, Ordering::Relaxed);
        if let Some(handle) = self.handle.take() {
            handle.join().ok();
        }
    }
}

/// A signal to a child we have not reaped yet always lands, so there is nothing
/// to check here: the pid cannot have been reused while we still owe it a wait.
fn signal(pid: libc::pid_t, sig: libc::c_int) {
    unsafe { libc::kill(pid, sig) };
}

/// Keep a failure's stderr and nothing else.
///
/// An empty file is not worth keeping either way, which is also how a command
/// that never spawned cleans up after itself: the file got created, the process
/// never wrote to it.
fn tidy(path: &Path, status: &Status) {
    let keep = status.failed() && std::fs::metadata(path).is_ok_and(|m| m.len() > 0);
    if !keep {
        std::fs::remove_file(path).ok();
    }
}

/// Local time as `20260813-142207`, for naming a directory after the run that
/// made it. Falls back to epoch seconds if the C library declines.
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

/// Start the command, handing back its pid and the moment it started. The
/// caller reaps it, which is what lets somebody else signal it in between.
fn spawn(cmd: &Cmd) -> anyhow::Result<(libc::pid_t, Instant)> {
    let mut proc = Command::new(&cmd.argv[0]);
    proc.args(&cmd.argv[1..]);
    proc.envs(&cmd.env);
    if let Some(dir) = &cmd.dir {
        proc.current_dir(dir);
    }
    proc.stdout(stdio(&cmd.stdout)?);
    proc.stderr(stdio(&cmd.stderr)?);

    let start = Instant::now();
    let child = proc
        .spawn()
        .with_context(|| format!("failed to spawn {}", cmd.argv[0]))?;

    Ok((child.id() as libc::pid_t, start))
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

/// Reap this specific child with `wait4` so its rusage comes back on its own.
/// `getrusage(RUSAGE_CHILDREN)` would fold in every sibling, which stops being
/// right the moment anything runs more than one process.
fn wait_timed(pid: libc::pid_t, start: Instant) -> anyhow::Result<Timing> {
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
