//! A pipeline end to end: a bare command, a batch, a serial step, two sinks.
//!
//! `cargo run -p pipeline --example basic`

use pipeline::{Cmd, OnError, Output, Pipeline, Progress, Step, Table};

fn main() -> anyhow::Result<()> {
    let dir = std::env::temp_dir().join("pipeline-basic");
    std::fs::remove_dir_all(&dir).ok();
    let data = dir.join("data.txt");
    let table = dir.join("runs.tbl");

    // four commands that each burn some cpu, run two at a time: the step's wall
    // clock comes out near double one command's, and the summed user time shows
    // what the parallelism bought
    let burns = (1..=4).map(|i| {
        Cmd::new(format!("burn-{i}"), "/bin/sh")
            .args(["-c", "seq 1 20000000 > /dev/null"])
            .label("job", i)
    });

    Pipeline::new()
        // a bare Cmd is a step of one
        .step(
            Cmd::new("make-data", "/bin/sh")
                .args(["-c", "seq 1 2000"])
                .stdout(Output::File(data.clone()))
                .tag("setup"),
        )
        .step(Step::batch("burn", 2, burns))
        .step(
            Step::serial(
                "count",
                [
                    Cmd::new("lines", "/usr/bin/wc")
                        .arg("-l")
                        .path(&data)
                        .stdout(Output::File(dir.join("lines")))
                        .label("mode", "lines"),
                    Cmd::new("bytes", "/usr/bin/wc")
                        .arg("-c")
                        .path(&data)
                        .stdout(Output::File(dir.join("bytes")))
                        .label("mode", "bytes"),
                    Cmd::new("missing-tool", "/no/such/tool"),
                ],
            )
            // one bad command should not cost us the rest of the matrix
            .on_error(OnError::Continue),
        )
        .sink(Progress::new())
        .sink(Table::new(&table))
        .run()?;

    println!("\n{}", table.display());
    print!("{}", std::fs::read_to_string(&table)?);

    Ok(())
}
