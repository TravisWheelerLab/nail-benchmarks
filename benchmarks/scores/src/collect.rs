//! Collecting a pipeline's shards in parallel and writing them in order.
//!
//! Every target lives in exactly one shard, so a shard's rows are a run of the
//! file that nothing outside it belongs in. Shards are collected in parallel
//! and written in shard order, which is what gets a sort of four billion rows
//! for the price of sorting each shard's own.
//!
//! The same bytes come out however many threads there were: a worker renders
//! its block into memory and the file takes them in shard order.
//!
//! Nothing here depends on what a row looks like. What a block holds is the
//! caller's, handed in as the renderer; what this owns is which shard is
//! collected when, how much memory that is allowed to cost, and the order the
//! blocks reach the file in.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex, mpsc};

use util::manifest;

use crate::shard::{Count, Scratch, Shard};

/// What a shard costs to hold while it is collected, beside the share of its
/// input that [`estimate`] allows for.
const OVERHEAD: u64 = 64 << 20;

/// How much of a shard's input it holds at the peak, in eighths.
//
// measured rather than reasoned: a synthetic shard of 4.0 GB
// peaked at 3.1 GB resident with one worker -- the hits and the
// domains it was read into, the block it was rendered to, and
// what the allocator held while each of those doubled
const PEAK: u64 = 6;

/// Which shards to collect, and what collecting them is allowed to cost.
pub struct Job<'a> {
    pub work: &'a Shard<'a>,
    pub shards: &'a [String],
    pub threads: usize,
    /// How many bytes the collectors may hold between them.
    pub mem: u64,
}

/// Collect every shard and write the blocks to `out` in shard order.
///
/// `block` renders one shard into the bytes that stand for it, on a worker
/// thread, and is the only part of this that varies by table.
pub fn blocks(
    job: Job<'_>,
    block: impl Fn(&Shard<'_>, &str, &mut Scratch) -> anyhow::Result<(Vec<u8>, Count)> + Sync,
    out: &mut impl Write,
) -> anyhow::Result<Count> {
    let Job {
        work,
        shards,
        threads,
        mem,
    } = job;

    let ticket = AtomicUsize::new(0);
    let stop = AtomicBool::new(false);
    let budget = Budget::new(mem);
    let (send, recv) = mpsc::channel::<Message>();

    let threads = threads.max(1).min(shards.len());

    std::thread::scope(|scope| -> anyhow::Result<Count> {
        for _ in 0..threads {
            let send = send.clone();
            let (ticket, stop, budget, block) = (&ticket, &stop, &budget, &block);

            scope.spawn(move || {
                let mut scratch = Scratch::default();

                loop {
                    if stop.load(Ordering::Relaxed) {
                        return;
                    }

                    let at = ticket.fetch_add(1, Ordering::Relaxed);
                    let Some(shard) = shards.get(at) else {
                        return;
                    };

                    let held = estimate(work, shard);
                    budget.take(at, held);

                    let message = match block(work, shard, &mut scratch) {
                        Ok((block, count)) => Message::Block {
                            at,
                            block,
                            count,
                            held,
                        },
                        Err(error) => {
                            budget.give(held);
                            stop.store(true, Ordering::Relaxed);
                            Message::Failed {
                                shard: shard.clone(),
                                error,
                            }
                        }
                    };

                    if send.send(message).is_err() {
                        return;
                    }
                }
            });
        }

        // the senders the workers hold are clones; this one would keep the
        // channel open after they are all done
        drop(send);

        let mut pending: BTreeMap<usize, (Vec<u8>, u64)> = BTreeMap::new();
        let mut next = 0usize;
        let mut total = Count::default();
        let mut failure: Option<(String, anyhow::Error)> = None;

        for message in recv {
            match message {
                Message::Failed { shard, error } => {
                    failure.get_or_insert((shard, error));

                    // nothing more will be written, so the blocks waiting for
                    // a shard that failed give their room back rather than
                    // holding the workers still reading
                    for (_, (_, held)) in std::mem::take(&mut pending) {
                        budget.give(held);
                    }
                }
                Message::Block {
                    at,
                    block,
                    count,
                    held,
                } => {
                    if failure.is_some() {
                        budget.give(held);
                        continue;
                    }

                    total.add(count);
                    pending.insert(at, (block, held));

                    while let Some((block, held)) = pending.remove(&next) {
                        out.write_all(&block)?;
                        next += 1;

                        budget.waiting_for(next);
                        budget.give(held);
                    }
                }
            }
        }

        if let Some((shard, error)) = failure {
            return Err(error.context(format!("shard {shard}")));
        }

        Ok(total)
    })
}

enum Message {
    Block {
        at: usize,
        block: Vec<u8>,
        count: Count,
        held: u64,
    },
    Failed {
        shard: String,
        error: anyhow::Error,
    },
}

/// What a shard will take to collect: its tables, the hits they become, and
/// the block they are rendered into.
fn estimate(work: &Shard<'_>, shard: &str) -> u64 {
    let seeds: u64 = work
        .lists
        .iter()
        .map(|list| {
            std::fs::metadata(manifest::seeds_path(work.results, list, shard))
                .map(|meta| meta.len())
                .unwrap_or(0)
        })
        .sum();

    let bytes: u64 = work
        .runs
        .iter()
        .filter(|column| column.shards.iter().any(|covered| covered == shard))
        .map(|column| {
            let table = manifest::table_path(work.results, &column.run.name, shard);
            let dom = manifest::dom_path(work.results, &column.run.name, shard);

            [table, dom]
                .iter()
                .filter_map(|path| std::fs::metadata(path).ok())
                .map(|meta| meta.len())
                .sum::<u64>()
        })
        .sum();

    (bytes + seeds) * PEAK / 8 + OVERHEAD
}

/// How much memory the collectors may hold between them.
struct Budget {
    limit: u64,
    state: Mutex<State>,
    room: Condvar,
}

struct State {
    used: u64,
    /// Which shard the file is waiting for.
    next: usize,
}

impl Budget {
    fn new(limit: u64) -> Budget {
        Budget {
            limit: limit.max(1),
            state: Mutex::new(State { used: 0, next: 0 }),
            room: Condvar::new(),
        }
    }

    fn take(&self, at: usize, want: u64) {
        let mut state = self.state.lock().expect("the budget outlives its holders");

        // three conditions, each admitting a shard that would
        // otherwise never run:
        //   used > 0        an oversized shard runs alone
        //                   rather than not at all
        //   +want <= limit  the ordinary case
        //   at == next      the shard the file is waiting for
        //
        // without the last the budget fills with blocks that are
        // collected and cannot be written: the file is waiting
        // for shard 4, the room is held by 5, 6 and 7, and the
        // worker that claimed 4 waits for room that only writing
        // 4 frees
        while state.used > 0 && state.used + want > self.limit && at != state.next {
            state = self
                .room
                .wait(state)
                .expect("the budget outlives its holders");
        }

        state.used += want;
    }

    fn give(&self, want: u64) {
        let mut state = self.state.lock().expect("the budget outlives its holders");
        state.used = state.used.saturating_sub(want);

        self.room.notify_all();
    }

    /// Say which shard the file is waiting for, which is the one that cannot
    /// be made to wait.
    fn waiting_for(&self, at: usize) {
        let mut state = self.state.lock().expect("the budget outlives its holders");
        state.next = at;

        self.room.notify_all();
    }
}

/// The machine's total memory, or eight gigabytes where it cannot be read.
pub fn ram() -> u64 {
    const FALLBACK: u64 = 8 << 30;

    #[cfg(target_os = "linux")]
    let total = std::fs::read_to_string("/proc/meminfo").ok().and_then(|text| {
        text.lines()
            .find_map(|line| line.strip_prefix("MemTotal:"))
            .and_then(|rest| rest.split_whitespace().next())
            .and_then(|kb| kb.parse::<u64>().ok())
            .map(|kb| kb * 1024)
    });

    #[cfg(target_os = "macos")]
    let total = std::process::Command::new("sysctl")
        .args(["-n", "hw.memsize"])
        .output()
        .ok()
        .and_then(|out| String::from_utf8(out.stdout).ok())
        .and_then(|text| text.trim().parse::<u64>().ok());

    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    let total: Option<u64> = None;

    total.unwrap_or(FALLBACK)
}

/// A path as the `#= cutoffs` line should carry it: what was given, made
/// absolute, so a table read elsewhere still names the file it was judged by.
pub fn absolute(path: &Path) -> PathBuf {
    match path.is_absolute() {
        true => path.to_path_buf(),
        false => std::env::current_dir()
            .map(|dir| dir.join(path))
            .unwrap_or_else(|_| path.to_path_buf()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_budget_admits_one_shard_bigger_than_all_of_it() {
        let budget = Budget::new(1 << 20);

        // nothing is held, so an oversized shard goes ahead rather than
        // waiting for room that will never come
        budget.take(0, 1 << 30);
        budget.give(1 << 30);
    }

    /// The shard the file is waiting for goes ahead even with the room full,
    /// since everything held is waiting on it.
    #[test]
    fn a_budget_never_holds_up_the_shard_being_written() {
        let budget = std::sync::Arc::new(Budget::new(1 << 20));
        budget.take(5, 1 << 20);
        budget.waiting_for(4);

        let waiting = budget.clone();
        let admitted = std::thread::spawn(move || waiting.take(4, 1 << 20));

        // join returns with the budget still full: nothing was released
        admitted.join().unwrap();
    }
}
