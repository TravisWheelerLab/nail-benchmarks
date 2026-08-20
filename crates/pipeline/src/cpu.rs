//! Which CPUs a command is allowed to run on.
//!
//! A command asks for a number of cores and never says which. [`Cores`] keeps
//! track of what is free, hands out that many, and takes them back when the
//! command is done; `taskset` does the pinning.

use std::sync::{Condvar, Mutex};

/// The cores a pipeline has to hand out, and which of them are in use.
///
/// The pool holds one logical CPU per physical core. Two logical CPUs on one
/// physical core are not two cores — they share the execution units, so a
/// command given both would get somewhere around 1.3 cores' worth of work done
/// while the table claimed it had 2. Only the first of each sibling group is
/// ever handed out and the rest sit idle, which is why a 32-CPU machine with
/// hyperthreading has 24 of these rather than 32.
#[derive(Debug, Default)]
pub(crate) struct Cores {
    pool: Vec<usize>,
    /// The cpus currently leased out, and something to wait on for one to come
    /// back. A command that cannot be placed yet sleeps here rather than
    /// spinning, which would burn a core to wait for a core.
    taken: Mutex<Vec<usize>>,
    freed: Condvar,
}

/// Cores held for as long as one command needs them.
///
/// Handing them back is [`Drop`]'s job rather than the caller's, so a command
/// that failed to spawn, timed out, or panicked releases its cores the same way
/// one that finished does.
pub(crate) struct Lease<'a> {
    cores: &'a Cores,
    cpus: Vec<usize>,
}

impl Cores {
    pub(crate) fn read() -> Cores {
        let mut pool = Vec::new();
        let mut spoken_for = Vec::new();

        for cpu in allowed() {
            if spoken_for.contains(&cpu) {
                continue;
            }
            spoken_for.extend(siblings(cpu));
            pool.push(cpu);
        }

        Cores {
            pool,
            taken: Mutex::new(Vec::new()),
            freed: Condvar::new(),
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.pool.len()
    }

    /// Take `size` cores, waiting for them if they are not free yet.
    ///
    /// `None` for a command that asked for none, and for one that gave up
    /// waiting because the run is stopping — neither is going to be pinned.
    /// A request larger than the whole machine would wait forever, which is why
    /// the pipeline refuses to build one.
    pub(crate) fn acquire(&self, size: usize, abandon: &dyn Fn() -> bool) -> Option<Lease<'_>> {
        if size == 0 || size > self.pool.len() {
            return None;
        }

        let mut taken = self.taken.lock().unwrap();
        loop {
            if abandon() {
                return None;
            }
            if let Some(lease) = self.grab(&mut taken, size) {
                return Some(lease);
            }
            taken = self.freed.wait(taken).unwrap();
        }
    }

    /// Take `size` cores if they are free right now, and never wait. For
    /// working out what a run would look like without running it.
    pub(crate) fn try_acquire(&self, size: usize) -> Option<Lease<'_>> {
        if size == 0 || size > self.pool.len() {
            return None;
        }
        self.grab(&mut self.taken.lock().unwrap(), size)
    }

    fn grab(&self, taken: &mut Vec<usize>, size: usize) -> Option<Lease<'_>> {
        // lowest first, so a run with the machine to itself places its commands
        // the same way every time
        let free: Vec<usize> = self
            .pool
            .iter()
            .filter(|cpu| !taken.contains(cpu))
            .take(size)
            .copied()
            .collect();

        if free.len() < size {
            return None;
        }

        taken.extend(&free);
        Some(Lease {
            cores: self,
            cpus: free,
        })
    }

    /// Wake anything waiting on cores that are never coming, so a stopping run
    /// does not leave workers parked.
    pub(crate) fn wake(&self) {
        self.freed.notify_all();
    }

    fn release(&self, cpus: &[usize]) {
        let mut taken = self.taken.lock().unwrap();
        taken.retain(|cpu| !cpus.contains(cpu));
        drop(taken);
        self.freed.notify_all();
    }
}

impl Lease<'_> {
    pub(crate) fn cpus(&self) -> &[usize] {
        &self.cpus
    }
}

impl Drop for Lease<'_> {
    fn drop(&mut self) {
        self.cores.release(&self.cpus);
    }
}

/// The wrapper that pins a command to `cpus`, or nothing at all if it was never
/// given any.
///
/// `taskset` execs into the program rather than supervising it, so the pid we
/// end up waiting on is the program's own and every number we measure is still
/// the program's. (`numactl` belongs here too, for a machine with more than one
/// node: it would want the node its cpus sit on, via `--membind`.)
pub(crate) fn wrapper(cpus: &[usize]) -> Vec<String> {
    if cpus.is_empty() {
        return Vec::new();
    }

    let list: Vec<String> = cpus.iter().map(|cpu| cpu.to_string()).collect();
    vec!["taskset".to_string(), "-c".to_string(), list.join(",")]
}

/// The CPUs this process may run on, which is not the same as every CPU the
/// machine has — a benchmark started under `taskset` gets a smaller set, and
/// handing out anything outside it would pin commands nowhere.
fn allowed() -> Vec<usize> {
    let mut set: libc::cpu_set_t = unsafe { std::mem::zeroed() };
    let rc = unsafe { libc::sched_getaffinity(0, size_of::<libc::cpu_set_t>(), &mut set) };
    if rc != 0 {
        return Vec::new();
    }

    (0..libc::CPU_SETSIZE as usize)
        .filter(|cpu| unsafe { libc::CPU_ISSET(*cpu, &set) })
        .collect()
}

/// Every logical CPU sharing a physical core with this one, itself included.
/// A CPU whose topology we cannot read is treated as a core of its own, which
/// is the reading that under-promises.
fn siblings(cpu: usize) -> Vec<usize> {
    let path = format!("/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list");
    match std::fs::read_to_string(path) {
        Ok(text) => parse_list(&text),
        Err(_) => vec![cpu],
    }
}

/// A kernel cpulist: `0-1`, `16`, `0-3,8-11`.
fn parse_list(text: &str) -> Vec<usize> {
    let mut out = Vec::new();

    for part in text.trim().split(',').filter(|p| !p.is_empty()) {
        match part.split_once('-') {
            Some((lo, hi)) => {
                if let (Ok(lo), Ok(hi)) = (lo.parse::<usize>(), hi.parse::<usize>()) {
                    out.extend(lo..=hi);
                }
            }
            None => out.extend(part.parse::<usize>()),
        }
    }
    out
}
