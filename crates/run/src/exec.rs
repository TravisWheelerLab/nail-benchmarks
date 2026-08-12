use std::fs::{File, OpenOptions};
use std::io;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::time::Instant;

use anyhow::{bail, Context};

/// Where a command's stdout or stderr goes.
///
/// Some tools write their results to stdout rather than to a named output file,
/// and some are invoked repeatedly with output appended to one table.
#[derive(Clone, Debug)]
pub enum Redirect {
    Null,
    Create(PathBuf),
    Append(PathBuf),
}

/// One process to spawn.
#[derive(Clone, Debug)]
pub struct Cmd {
    pub argv: Vec<String>,
    pub stdout: Redirect,
    pub stderr: Redirect,
}

impl Cmd {
    pub fn new(program: impl AsRef<std::path::Path>) -> Self {
        Cmd {
            argv: vec![program.as_ref().display().to_string()],
            stdout: Redirect::Null,
            stderr: Redirect::Null,
        }
    }

    pub fn arg(mut self, arg: impl Into<String>) -> Self {
        self.argv.push(arg.into());
        self
    }

    /// A path argument, spelled out so callers do not repeat `.display()`.
    pub fn path(self, path: impl AsRef<std::path::Path>) -> Self {
        self.arg(path.as_ref().display().to_string())
    }

    pub fn args<I, S>(mut self, args: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.argv.extend(args.into_iter().map(Into::into));
        self
    }

    pub fn stdout_to(mut self, path: impl Into<PathBuf>) -> Self {
        self.stdout = Redirect::Create(path.into());
        self
    }

    pub fn stdout_append(mut self, path: impl Into<PathBuf>) -> Self {
        self.stdout = Redirect::Append(path.into());
        self
    }

    /// Capture stderr rather than letting it reach the console. Tools are
    /// chatty, but their output is the only diagnostic when a run fails.
    pub fn stderr_to(mut self, path: impl Into<PathBuf>) -> Self {
        self.stderr = Redirect::Append(path.into());
        self
    }
}

/// The full argv as a single line, for the `cmd` column of runs.tbl.
pub fn render(cmd: &Cmd, numa: Option<&Numa>) -> String {
    let mut parts: Vec<String> = Vec::new();
    if let Some(numa) = numa {
        parts.extend(numa.prefix());
    }
    parts.extend(cmd.argv.iter().cloned());
    parts.join(" ")
}

fn command(cmd: &Cmd, numa: Option<&Numa>) -> anyhow::Result<Command> {
    let argv = match numa {
        Some(numa) => {
            let mut v = numa.prefix();
            v.extend(cmd.argv.iter().cloned());
            v
        }
        None => cmd.argv.clone(),
    };

    let mut out = Command::new(&argv[0]);
    out.args(&argv[1..]);
    out.stdout(stdio(&cmd.stdout)?);
    out.stderr(stdio(&cmd.stderr)?);
    Ok(out)
}

fn stdio(out: &Redirect) -> anyhow::Result<Stdio> {
    Ok(match out {
        Redirect::Null => Stdio::null(),
        Redirect::Create(path) => {
            if let Some(dir) = path.parent() {
                std::fs::create_dir_all(dir)?;
            }
            let file = File::create(path)
                .with_context(|| format!("failed to create {}", path.display()))?;
            Stdio::from(file)
        }
        Redirect::Append(path) => {
            if let Some(dir) = path.parent() {
                std::fs::create_dir_all(dir)?;
            }
            let file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)
                .with_context(|| format!("failed to open {} for append", path.display()))?;
            Stdio::from(file)
        }
    })
}

/// Resource usage for one child process, taken from `wait4`.
///
/// The numbers come straight out of the kernel rather than being formatted to
/// text and parsed back.
#[derive(Clone, Copy, Debug)]
pub struct Timing {
    pub wall_s: f64,
    pub user_s: f64,
    pub sys_s: f64,
    pub max_rss_kb: i64,
    pub exit: i32,
}

impl Timing {
    /// Aggregate several children into one logical run: wall clock is the batch
    /// span, CPU time sums, peak RSS is the high-water mark, and the exit code
    /// is the first nonzero one.
    pub fn combine(parts: &[Timing], wall_s: f64) -> Timing {
        Timing {
            wall_s,
            user_s: parts.iter().map(|p| p.user_s).sum(),
            sys_s: parts.iter().map(|p| p.sys_s).sum(),
            max_rss_kb: parts.iter().map(|p| p.max_rss_kb).max().unwrap_or(0),
            exit: parts.iter().map(|p| p.exit).find(|c| *c != 0).unwrap_or(0),
        }
    }
}

/// CPU pinning via `numactl`. Only constructed when a config actually asks for
/// a node, so a non-NUMA machine never invokes `numactl` at all.
#[derive(Clone, Debug)]
pub struct Numa {
    node: usize,
    cpus: Vec<usize>,
}

impl Numa {
    /// Pin to the first `threads` CPUs of `node`.
    pub fn new(node: usize, threads: usize) -> anyhow::Result<Self> {
        let out = Command::new("numactl")
            .arg("--hardware")
            .output()
            .context("failed to run `numactl --hardware`; is numactl installed?")?;

        if !out.status.success() {
            bail!("`numactl --hardware` failed with status {}", out.status);
        }

        let text = String::from_utf8_lossy(&out.stdout);
        let line = text
            .lines()
            .find(|l| {
                let f: Vec<&str> = l.split_whitespace().collect();
                f.len() > 3 && f[0] == "node" && f[1] == node.to_string() && f[2] == "cpus:"
            })
            .with_context(|| format!("numa node {node} not present in `numactl --hardware`"))?;

        let cpus: Vec<usize> = line
            .split_whitespace()
            .skip(3)
            .filter_map(|c| c.parse().ok())
            .take(threads)
            .collect();

        if cpus.len() < threads {
            bail!(
                "numa node {node} exposes {} cpus, fewer than the {threads} requested",
                cpus.len()
            );
        }

        Ok(Numa { node, cpus })
    }

    fn prefix(&self) -> Vec<String> {
        let list = self
            .cpus
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join(",");

        vec![
            "numactl".to_string(),
            format!("--physcpubind={list}"),
            format!("--membind={}", self.node),
        ]
    }
}

/// Run one command to completion, returning its resource usage.
pub fn run(cmd: &Cmd, numa: Option<&Numa>) -> anyhow::Result<Timing> {
    let mut command = command(cmd, numa)?;
    let start = Instant::now();

    let child = command
        .spawn()
        .with_context(|| format!("failed to spawn {}", cmd.argv[0]))?;

    wait_timed(child.id() as libc::pid_t, start)
}

/// Run several commands at once and aggregate them into a single Timing.
pub fn run_together(cmds: &[Cmd], numa: Option<&Numa>) -> anyhow::Result<Timing> {
    let start = Instant::now();

    let mut pids = Vec::with_capacity(cmds.len());
    for cmd in cmds {
        let mut command = command(cmd, numa)?;
        let child = command
            .spawn()
            .with_context(|| format!("failed to spawn {}", cmd.argv[0]))?;
        pids.push(child.id() as libc::pid_t);
    }

    let mut parts = Vec::with_capacity(pids.len());
    for pid in pids {
        parts.push(wait_timed(pid, start)?);
    }

    Ok(Timing::combine(&parts, start.elapsed().as_secs_f64()))
}

/// Run several commands in turn and aggregate them into a single Timing.
pub fn run_each(cmds: &[Cmd], numa: Option<&Numa>) -> anyhow::Result<Timing> {
    let start = Instant::now();

    let mut parts = Vec::with_capacity(cmds.len());
    for cmd in cmds {
        parts.push(run(cmd, numa)?);
    }

    Ok(Timing::combine(&parts, start.elapsed().as_secs_f64()))
}

/// Run a command that is not part of a measurement and fail loudly if it does
/// not succeed.
pub fn check(cmd: &Cmd, numa: Option<&Numa>, what: &str) -> anyhow::Result<()> {
    let timing = run(cmd, numa)?;
    if timing.exit != 0 {
        bail!(
            "{what} failed with exit code {}: {}",
            timing.exit,
            render(cmd, numa)
        );
    }
    Ok(())
}

/// Reap a specific child via wait4 so we get its rusage in isolation.
/// getrusage(RUSAGE_CHILDREN) would accumulate across siblings, which is wrong
/// as soon as we run commands concurrently.
fn wait_timed(pid: libc::pid_t, start: Instant) -> anyhow::Result<Timing> {
    let mut status: libc::c_int = 0;
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };

    // std::process::Child has no Drop that reaps, so taking the child with
    // wait4 here does not race with anything in std.
    let rc = unsafe { libc::wait4(pid, &mut status, 0, &mut usage) };
    if rc < 0 {
        return Err(io::Error::last_os_error()).context("wait4 failed");
    }

    let wall_s = start.elapsed().as_secs_f64();

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn timing_captures_cpu_and_exit_code() {
        let cmd = Cmd::new("/bin/sh").arg("-c").arg("exit 3");
        let timing = run(&cmd, None).unwrap();
        assert_eq!(timing.exit, 3);
        assert!(timing.wall_s >= 0.0);
        assert!(timing.max_rss_kb > 0, "expected nonzero peak rss");
    }

    #[test]
    fn signalled_child_reports_128_plus_signal() {
        let cmd = Cmd::new("/bin/sh").arg("-c").arg("kill -TERM $$");
        let timing = run(&cmd, None).unwrap();
        assert_eq!(timing.exit, 128 + libc::SIGTERM);
    }

    #[test]
    fn stdout_can_be_redirected_to_a_file() {
        let dir = std::env::temp_dir().join(format!("bm-exec-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let out = dir.join("out.txt");

        let cmd = Cmd::new("/bin/sh")
            .arg("-c")
            .arg("echo hello")
            .stdout_to(&out);
        let timing = run(&cmd, None).unwrap();

        assert_eq!(timing.exit, 0);
        assert_eq!(std::fs::read_to_string(&out).unwrap().trim(), "hello");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn concurrent_commands_sum_cpu_but_not_wall_time() {
        let cmds: Vec<Cmd> = (0..3)
            .map(|_| {
                Cmd::new("/bin/sh")
                    .arg("-c")
                    .arg("i=0; while [ $i -lt 200000 ]; do i=$((i+1)); done")
            })
            .collect();

        let timing = run_together(&cmds, None).unwrap();
        assert_eq!(timing.exit, 0);
        // three busy children should burn more cpu than the batch took in wall
        // clock, which is the whole reason we sum user time rather than wall
        assert!(
            timing.user_s + timing.sys_s > timing.wall_s,
            "cpu {} should exceed wall {}",
            timing.user_s + timing.sys_s,
            timing.wall_s
        );
    }
}
