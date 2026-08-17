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

const WATCHDOG_TICK: Duration = Duration::from_millis(100);
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
}


pub fn execute(cmd: &mut Cmd) {
    let live = Arc::new(Live::default());
    // watch if there's a timeout
    let _watchdog = cmd.timeout.is_some().then(|| Live::watch(&live));
    live.execute(cmd);
}

#[derive(Default)]
pub(crate) struct Live {
    running: Mutex<Vec<Entry>>,
    stopping: AtomicBool,
    done: AtomicBool,
}

struct Entry {
    pid: libc::pid_t,
    deadline: Option<Instant>,
    phase: Phase,
    expired: bool,
}

enum Phase {
    Running,
    Term(Instant),
    Kill,
}

impl Live {
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

    pub(crate) fn stop(&self) {
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

    pub(crate) fn watch(live: &Arc<Live>) -> Watchdog {
        let live = Arc::clone(live);
        let handle = std::thread::spawn({
            let live = Arc::clone(&live);
            move || {
                while !live.done.load(Ordering::Relaxed) {
                    live.sweep();
                    std::thread::sleep(WATCHDOG_TICK);
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

    fn remove(&self, pid: libc::pid_t) -> bool {
        let mut running = self.running.lock().unwrap();
        match running.iter().position(|e| e.pid == pid) {
            Some(i) => running.swap_remove(i).expired,
            None => false,
        }
    }

    fn sweep(&self) {
        let now = Instant::now();
        let mut running = self.running.lock().unwrap();
        for entry in running.iter_mut() {
            match entry.phase {
                Phase::Running if entry.deadline.is_some_and(|at| now >= at) => {
                    entry.expired = true;
                    entry.term(now);
                }
                Phase::Term(at) if now.duration_since(at) >= SIGTERM_GRACE => {
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

fn signal(pid: libc::pid_t, sig: libc::c_int) {
    unsafe { libc::kill(pid, sig) };
}

fn tidy(path: &Path, status: &Status) {
    let keep = status.failed() && std::fs::metadata(path).is_ok_and(|m| m.len() > 0);
    if !keep {
        std::fs::remove_file(path).ok();
    }
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

fn spawn(cmd: &Cmd) -> anyhow::Result<(libc::pid_t, Instant)> {
    let mut proc = Command::new(&cmd.program);
    proc.args(cmd.args());
    proc.envs(&cmd.env);
    if let Some(dir) = &cmd.dir {
        proc.current_dir(dir);
    }
    proc.stdout(stdio(&cmd.stdout)?);
    proc.stderr(stdio(&cmd.stderr)?);

    let start = Instant::now();
    let child = proc
        .spawn()
        .with_context(|| format!("failed to spawn {}", cmd.program.display()))?;

    Ok((child.id() as libc::pid_t, start))
}

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
