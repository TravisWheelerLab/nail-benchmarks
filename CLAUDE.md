# nail-benchmarks

Benchmarks for [nail](https://github.com/travisWheelerLab/nail), a profile HMM
search tool, against the tools it is compared to: HMMER, MMseqs2, BLAST, LAST
and DIAMOND.

A Cargo workspace at the repo root. Binaries are the interface. There is one
per benchmark, and each owns its own inputs, runs and analyses.

Each binary has a shim beside its `Cargo.toml`, named after it, so the
`mgy build` and `pid run` spellings used throughout this file are what you
type:

```
benchmarks/mgy/mgy recall --shards 2
benchmarks/pid/pid parse recall
benchmarks/long-seqs/long-seqs parse cells
```

All three are symlinks to `benchmarks/shim`; the crate it builds is the name
on the link. It builds that crate and runs it only if the build succeeded, so
a broken build never falls through to the last binary that worked. Nothing can
ask cargo whether it would rebuild without letting it rebuild, so the shim
does not try: it runs `cargo build` every time, which costs about a twentieth
of a second when there is nothing to do. A build still going after that gets a
spinner on a terminal, and nothing at all into a pipe or a log file.

On macOS the shim runs a copy of the binary kept under `target/shim/`,
refreshed only when cargo produced different bytes. cargo replaces
`target/release/<name>` on every build, no-op builds included, and macOS
spends a fifth of a second checking the signature of an executable whose
inode is new to it. Without the copy the shim would cost a quarter of a
second rather than a twentieth.

## Layout

```
Makefile               downloads data, builds tools, nothing else
data/                  what the Makefile downloaded
tools/bin/             what the Makefile built
benchmarks/shim        the build-and-run shim every benchmark links to
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
- `tabl` writes the padded, `#`-headed tables. It is not published: the
  workspace takes it as a path dependency on a sibling checkout, `../tabl`.
- `feisty` sits in `[workspace.dependencies]` and no member uses it.

## The shape a benchmark has

An input set under `inputs/`, one directory per pipeline under `outputs/`, and
the scratch off to one side:

```
inputs/<set>/                    what a run reads
outputs/<pipeline>/
├── manifest.tbl                 every command, its wall clock, its exit code
├── ledger.tbl                   one row per run per shard, and what it cost
└── results/<run>.<shard>.tbl    one hit table per run, per target shard
tmp/<pipeline>/                  scratch, and nothing worth keeping
```

The scratch sits at the crate root rather than inside the pipeline's output
directory, one subdirectory per pipeline: `mgy/tmp/recall`,
`mgy/tmp/cloud-search`, `mgy/tmp/cutoffs/<calibration>/<stage>`, and
`mgy/tmp/build` for what `build` writes on the way. So `outputs/<pipeline>/`
holds the record and only the record, and the whole of `tmp/` can go at any
time without touching it. The six search pipelines still take `--tmp` to put
their own scratch somewhere else, a scratch disk being the usual reason;
`build` and `cutoffs` have no such flag and never had one.

`<benchmark> clean` removes `inputs/`, `outputs/` and `tmp/`, after printing
how many files and how many bytes are about to go and waiting for a `y`. The
two expensive trees outside those names, mgy's `cutoffs/` and pid's
`profmark/`, go only with `--all`. long-seqs' `inputs/` are symlinks to what
is checked in under `data/long-seqs/`, so its clean leaves them where they
are.

`ledger.tbl` is what the analyses read, and `parse` reads a pipeline's shape
out of it rather than out of the filenames. That is what keeps the analyses
free of any one pipeline: a row's `name` becomes a column, its `tool` says how
to read that column's table, its `shard` names the target it searched, and
every other column is a setting that tells one run from another. Adding a run
to a sweep adds a column without touching the reader.

A pipeline writes its own ledger when it finishes, distilled out of the
`manifest.tbl` `michi` left beside it. That is where a batched step folds to
its longest command and a command that failed is left out, once rather than in
every analysis. A pipeline also clears the ledger before it starts, so a run
that dies partway leaves none. The analyses then fall back to distilling the
manifest of the run that just failed and say what did not finish, rather than
reading the ledger of the run before it and reporting numbers for results that
have since been overwritten.

Results that came back from a cluster have no manifest and cannot have one,
since nothing here ran those commands and an exit code or an argv would have to
be invented. They arrive with `.time` files beside the tables, and `mgy import`
reads those. It is the only command that writes a ledger by hand.

## benchmarks/util

- `manifest` reads back the table `michi`'s sink wrote, and builds the
  `results/` paths from it.
- `ledger` holds the shape the analyses read: one row per run per shard, with
  each run's seconds already totalled. It distills a manifest into that shape,
  and reads and writes `ledger.tbl`.
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

`mgy build` writes `sizes.tbl` beside the target shards as it deals them:
counting a thousand shards afterwards is the whole deal read again, and the
count is only a metadata line. `mgy build sizes` writes one for a set dealt
before that.

`mgy parse scores` reads recall into `scores.tbl`: one row per query/target
pair, one score column per tool, and a `pass` string holding one character per
run. It collects each shard on its own and writes it as a block in shard order,
so the file comes out the same bytes however many threads it ran; `--threads`
and `--mem` size that. `summary` is one streaming pass over the result.

Each of the other two pipelines gets its own table, its own `tabl` schema and
its own reader when it is built; the shard-parallel collector underneath is
shared. Until hit-loss has one, `funnel` reads the older shape and nothing
writes it.

Two things sit outside the shape above. `mgy cutoffs` is a calibration rather
than a benchmark: five stages that reverse the targets, recruit decoys per
family, search each family against its own decoys forward and reversed, and
learn the per-family score cutoffs every hit is then held against. It keeps its
own directory tree, and produces `data/mgy-cutoffs.tbl`, which is committed and
promoted by hand. `mgy import` writes a `ledger.tbl` for result tables produced
elsewhere, out of the `.time` files that came back with them, which turns a
search run on a cluster into an ordinary pipeline directory.
`benchmarks/mgy/scripts/rename-old-results.sh` renames the older harness's
files into the names it expects.

Two dev tools sit beside all of that and belong to no pipeline.
`src/bin/synth_results.rs` writes a recall directory of the right shape and
size without a search behind it, so `parse scores` can be timed against one the
size of a real run; it is the crate's second binary, and the shim does not run
it. `scripts/compare-scores.py` reduces a `scores.tbl` and one in the older
shape to the same sets of pairs, tool scores and pass flags, and says where
they differ.

## benchmarks/pid

Recall as a function of the percent identity between a query and its target.
`pid build` assembles the benchmark from a profmark split, holding both axes
and the truth table together: `benchmark.tbl` records which pair is which and
at what identity, and belongs to neither side. There is one benchmark, under
`inputs/`, and `build` refuses to overwrite it. The profmark split itself sits
at the crate root, since it is expensive, depends only on Pfam and the split
parameters, and a rebuild draws from the same one.

`pid run` searches every tool against the benchmark, `parse` turns the results
into the tables the plot scripts read, and `plot` draws them.

## benchmarks/long-seqs

How each tool's runtime scales as sequences get longer. Six paired
query/target files, where `run` searches each query against its pair and
`parse` turns that into the plot scripts' tables. Its inputs are small and
checked in under `data/long-seqs/`, so only its `outputs/` is ignored.

## Testing behaviour: work in your own copy

This working tree belongs to whoever is at the keyboard. A pipeline writes into
`benchmarks/*/outputs/`, a `build` subcommand rewrites `inputs/`, and when two
people write there at once neither can tell which results are theirs. So run
nothing here. Copy the source into `tmp-claude/sandbox/`, link the expensive
directories, and work in the copy.

```bash
ROOT=$(git rev-parse --show-toplevel)
SB=$ROOT/tmp-claude/sandbox

rsync -a --delete \
  --exclude '.git/' --exclude 'target/' --exclude 'tmp-claude/' \
  --exclude '/data' --exclude '/tools' --exclude '/benchmarks/pid/profmark' \
  --exclude 'outputs/' --exclude 'tmp/' \
  --exclude '/benchmarks/mgy/inputs/' --exclude '/benchmarks/pid/inputs/' \
  "$ROOT/" "$SB/"

ln -sfn "$ROOT/data" "$SB/data"
ln -sfn "$ROOT/tools" "$SB/tools"
ln -sfn "$ROOT/benchmarks/pid/profmark" "$SB/benchmarks/pid/profmark"

# the workspace names tabl `../tabl`, which from the copy is this
ln -sfn "$(dirname "$ROOT")/tabl" "$ROOT/tmp-claude/tabl"
```

Run that before testing anything, and again after every edit: a copy goes
stale, and a result from stale source is worth nothing. `--delete` is what
keeps it current, and it drops what was deleted from the source without
touching the sandbox's own `inputs/`, `outputs/`, `target/` or links, since
rsync leaves excluded paths on the receiving side alone. The three link paths are excluded
without a trailing slash on purpose: a pattern ending in `/` matches only a
directory, and on the receiving side these are symlinks, so `--delete` removes
them. It copies uncommitted edits, which is the point: what wants testing is
usually not committed yet.

Then run inside `$SB`, through its own shims, which build there and run what
they built. Every path these crates resolve comes from their own
`CARGO_MANIFEST_DIR` (`util/src/tools.rs:18`, `mgy/src/main.rs:84`, and the
same in pid and long-seqs), so a build in the sandbox reads the sandbox's
`data/` and `tools/` links and writes the sandbox's `outputs/`. Re-syncing an
existing copy is near instant, and the build from cold takes about fifteen
seconds.

The links are the two directories worth 4.3G between them, plus the profmark
split, which is expensive for the reasons the pid section gives. The sandbox
reads all three; nothing in it should write them, so do not run `make data` or
`make tools` from the copy. The fourth is source rather than data: `tabl` is
built from wherever that link points, so an edit to it reaches the sandbox
without a sync.

mgy's and pid's `inputs/` are excluded for the same reason `outputs/` is: a
build writes them and nothing commits them, so the copy builds its own. Without
the exclude, `--delete` removes the set a `mgy build` in the sandbox just
wrote, since the real tree has nothing there to match it. long-seqs' `inputs/`
are two checked-in symlinks into `data/`, so they are copied like any other
tracked file.

Edit in the real tree and re-sync, never in the sandbox: an edit in the copy is
gone at the next sync. To check an analysis against a run that finished
elsewhere, copy that run's `outputs/<pipeline>/` into the sandbox and parse it
there.

## What is tracked

The downloads under `data/` and the builds under `tools/bin/` are not, and
neither are the generated input sets or any pipeline's output. Two exceptions
are committed on purpose: `data/long-seqs/`, because nothing fetches it and
ignoring it would empty a fresh clone, and `data/mgy-cutoffs.tbl`, because
learning it is a calibration run rather than a download.

Plotting is Python and matplotlib, under each benchmark's `scripts/`.
