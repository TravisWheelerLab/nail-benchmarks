# nail-benchmarks

Benchmarks for [nail](https://github.com/travisWheelerLab/nail), a profile HMM
search tool, against the tools it is compared to: HMMER, MMseqs2, BLAST, LAST
and DIAMOND.

A Cargo workspace at the repo root. Binaries are the interface. There is one
per benchmark, and each owns its own inputs, runs and analyses.

## Layout

```
Makefile               downloads data, builds tools, nothing else
data/                  what the Makefile downloaded
tools/bin/             what the Makefile built
benchmarks/util/       what the benchmarks share
benchmarks/mgy/        Pfam against MGnify
benchmarks/pid/        recall against percent identity, over a profmark split
benchmarks/long-seqs/  how the tools scale with sequence length
```

No benchmark looks on `PATH`. `util::tools` holds the path to every binary and
every download, and a benchmark reads it rather than guessing, so a run uses
whichever `hmmsearch` was built for this repo.

## External crates

- `michi` runs a pipeline of commands, times each one, and writes what it did to
  `manifest.tbl`.
- `libsail` reads and writes the formats: FASTA, Stockholm, p7hmm, and the hit
  tables nail, HMMER, MMseqs2 and BLAST produce.
- `feisty` sits in `[workspace.dependencies]` and no member uses it.

## The shape a benchmark has

An input set under `inputs/`, and one directory per pipeline under `outputs/`:

```
inputs/<set>/                    what a run reads
outputs/<pipeline>/
├── manifest.tbl                 every command, its wall clock, its exit code
├── results/<run>.<shard>.tbl    one hit table per run, per target shard
└── tmp/                         scratch
```

`manifest.tbl` is the record of what ran, and `parse` reads it instead of the
filenames. That is what keeps the analyses free of any one pipeline: a row
carrying a `name` field becomes a column, its `tool` field says how to read its
table, and its `shard` field names the target it searched. Adding a run to a
sweep adds a column without touching the reader.

## benchmarks/util

- `manifest` reads back the table `michi`'s sink wrote, and builds the
  `results/` paths from it.
- `tbl` writes the padded, `#`-headed table every analysis produces.
- `tools` holds where the binaries and the downloads are.
- `split` cuts a query set into balanced parts for a batch of jobs.
- `nail` reads the one column of nail's table that `libsail`'s layout does not
  carry.
- `time` reads the `.time` file of a run made outside the harness, in any of
  the six formats a `time` command might have written.

## benchmarks/mgy

Pfam profiles against MGnify metagenomic sequences, and the largest of the
three. `mgy build` cuts the sources into an input set: `fixed` is one query set
against target shards of equal size, `ladder` is nested rungs on both axes,
each a prefix of the one above.

Four pipelines. `recall` searches every shard while sweeping nail's and
MMseqs2's prefilter sensitivity. `cloud-search` seeds once, then searches every
`(A, B)` pruning cell off those same seeds, so the pruning parameters are the
only thing moving. `hit-loss` follows where the hits HMMER finds get lost.
`search-size` times every tool over every rung of both ladders.

`mgy parse scores` reads any of them into `scores.tbl`, which holds one row per
query/target pair and one column per run. `summary` and `funnel` are groupings
over that.

Two things sit outside the shape above. `mgy cutoffs` is a calibration rather
than a benchmark: five stages that reverse the targets, recruit decoys per
family, search each family against its own decoys forward and reversed, and
learn the per-family score cutoffs every hit is then held against. It keeps its
own directory tree, and produces `data/mgy-cutoffs.tbl`, which is committed and
promoted by hand. `mgy install` writes a `manifest.tbl` for result tables
produced elsewhere, which turns a search run on a cluster into an ordinary
pipeline directory. `benchmarks/mgy/scripts/rename-old-results.sh` renames the
older harness's files into the names it expects.

## benchmarks/pid

Recall as a function of the percent identity between a query and its target.
`pid build` assembles an input set from a profmark split, holding both axes and
the truth table together: `benchmark.tbl` records which pair is which and at
what identity, and belongs to neither side. Input sets are named by size, and
the profmark split itself sits at the crate root, since it is expensive,
depends only on Pfam and the split parameters, and every size is drawn from it.

`pid run` searches every tool against the benchmark, `parse` turns the results
into the tables the plot scripts read, and `plot` draws them.

## benchmarks/long-seqs

How each tool's runtime scales as sequences get longer. Six paired
query/target files, where `run` searches each query against its pair and
`parse` turns that into the plot scripts' tables. Its inputs are small and
checked in under `data/long-seqs/`, so only its `outputs/` is ignored.

## What is tracked

The downloads under `data/` and the builds under `tools/bin/` are not, and
neither are the generated input sets or any pipeline's output. Two exceptions
are committed on purpose: `data/long-seqs/`, because nothing fetches it and
ignoring it would empty a fresh clone, and `data/mgy-cutoffs.tbl`, because
learning it is a calibration run rather than a download.

Plotting is Python and matplotlib, under each benchmark's `scripts/`.
