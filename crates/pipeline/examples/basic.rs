//! A pipeline end to end: a bare command, a batch, a serial step, a deadline, a
//! step that gives up, and two sinks.
//!
//! `cargo run -p pipeline --example basic`

use std::time::Duration;

use pipeline::{Cmd, OnError, PipelineBuilder, Progress, Step, Table};

fn main() -> anyhow::Result<()> {
    let dir = std::env::temp_dir().join("pipeline-basic");
    std::fs::remove_dir_all(&dir).ok();
    let data = dir.join("data.txt");
    let table = dir.join("runs.tbl");

    // four commands that each burn some cpu, run two at a time: the step's wall
    // clock comes out near double one command's, and the summed user time shows
    // what the parallelism bought
    let burns = (1..=4).map(|i| {
        Cmd::new("/bin/sh")
            .name(format!("burn-{i}"))
            .arg("-c")
            // arg takes anything printable, so the count needs no to_string
            .arg(format_args!("seq 1 {} > /dev/null", i * 5_000_000))
            .field("job", i)
    });

    // three that would take half a minute each, alongside one that fails at
    // once. the step is set to stop, so the sleeps get killed rather than run to
    // the end for nothing
    let sleeps = (1..=2).map(|i| {
        Cmd::new("/bin/sh")
            .name(format!("sleep-{i}"))
            .args(["-c", "sleep 30"])
    });

    PipelineBuilder::new()
        // a bare Cmd is a step of one
        .step(
            Cmd::new("/bin/sh")
                .name("make-data")
                .args(["-c", "seq 1 2000"])
                .stdout_to(&data)
                .tag("setup"),
        )
        .step(Step::batch(2, burns).name("burn"))
        .step(
            Step::serial([
                // dir means the argument can be relative, and the variable rides
                // along in the pasteable line
                Cmd::new("/usr/bin/wc")
                    .name("lines")
                    .arg("-l")
                    .path("data.txt")
                    .dir(&dir)
                    .env("LC_ALL", "C")
                    .stdout_to(dir.join("lines"))
                    .field("mode", "lines"),
                Cmd::new("/usr/bin/wc")
                    .name("bytes")
                    .arg("-c")
                    .path("data.txt")
                    .dir(&dir)
                    .env("LC_ALL", "C")
                    .stdout_to(dir.join("bytes"))
                    .field("mode", "bytes"),
                // no name: it goes by its number alone
                Cmd::new("/no/such/tool"),
                // fails with something to say, so its stderr is kept
                Cmd::new("/bin/sh")
                    .name("boom")
                    .args(["-c", "echo 'no such database' >&2; exit 3"]),
            ])
            .name("count")
            // one bad command should not cost us the rest of the matrix
            .on_error(OnError::Continue),
        )
        .step(
            Step::serial([Cmd::new("/bin/sh")
                .name("slow")
                .args(["-c", "sleep 30"])
                .timeout(Duration::from_millis(300))])
            .name("limit")
            .on_error(OnError::Continue),
        )
        .step(
            Step::batch(
                3,
                std::iter::once(
                    Cmd::new("/bin/sh")
                        .name("early")
                        .args(["-c", "exit 1"])
                        .field("job", 0),
                )
                .chain(sleeps),
            )
            .name("race"),
        )
        // never gets a turn, because the step above stops the run
        .step(Cmd::new("/bin/echo").name("after").arg("unreachable"))
        .sink(Progress::new())
        .sink(Table::new(&table))
        .build()
        .run()?;

    println!("\n{}", table.display());
    print!("{}", std::fs::read_to_string(&table)?);

    Ok(())
}
