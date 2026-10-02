# nail-benchmarks

Benchmarks for [nail](https://github.com/travisWheelerLab/nail), a profile HMM
search tool, against the tools it is compared to: HMMER, MMseqs2, BLAST, LAST
and DIAMOND.

A Cargo workspace at the repo root. Binaries are the interface. There is one
per benchmark, and what each of them builds, runs and works out goes under
`sets/`.

Work and facts live in foam, not here: `foam ready` for what is open, `foam
memories` for what earlier sessions measured. This file holds the rules.

Each binary has a shim beside its `Cargo.toml`, named after it, so the
`recall run` and `pid parse` spellings used throughout this file are what you
type:

```
benchmarks/recall/recall run --in real --threads 32
benchmarks/build-set/build-set --in recall-toy --threads 8
benchmarks/pid/pid parse recall
```

Every one of them is a symlink to `benchmarks/shim`; the crate it builds is
the name on the link. It builds that crate and runs it only if the build
succeeded, so a broken build never falls through to the last binary that
worked. It runs `cargo build` every time, since nothing can ask cargo whether
it would rebuild without letting it rebuild. On macOS it runs a copy kept
under `target/shim/`, refreshed only when cargo produced different bytes,
because macOS re-checks the signature of any executable whose inode is new.

## Working with Jack

- Replies are a handful of lines: the result and the decision he owes. Keep
  evidence, measurements and the reasoning behind a recommendation for when he
  asks. A caveat that changes what he would decide still goes in, in a
  sentence.
- Answer a question in his own terms, yes or no or the fact first. Never open
  by telling him his question rests on a mistake or by restating it as
  something narrower.
- "Be careful" and "don't break my library" mean build the complete, general
  feature and verify it. They never mean scope it down to what the current
  consumer needs. Where a case is impossible, say so as a property of the
  format and offer a fix.
- Never put `rm -rf`, `rm -r`, `find -delete` or `rsync --delete` on a path
  assembled from shell variables. Name the literal path and look at it first.
  For experiments, write each run into a fresh, uniquely named directory
  rather than clearing one. Anything over a line or two goes in a script under
  `tmp-claude/` with `set -euo pipefail` and quoted variables, shown before it
  runs. Never write or delete inside `sets/`.
- Dropping a data point, or changing what a curve is measured against, is a
  decision for Jack up front, never a footnote in a deviations list.
- Delegate surveys across the crates to sub-agents and keep only the
  conclusion in the main context.

## Layout

```
Makefile                  downloads data, builds tools, nothing else
data/                     what the Makefile downloaded
tools/bin/                what the Makefile built
sets/                     what a build, a run or an analysis produced
sets-toy/                 the same, for the toy of each
figures/                  the pdfs, flat
figures-toy/              the same, for the toys
reports/                  hand-assembled results packages, untracked
benchmarks/shim           the build-and-run shim every binary links to

benchmarks/util/          set, search, paths, manifest, ledger, tools,
                          split, cut, clean, domains

benchmarks/build-set/     cuts sources into a set: fixed | reversed | cross |
                          pairs | profmark

benchmarks/recall/        recall against prefilter sensitivity   [fixed]
benchmarks/cloud-search/  the (A, B) pruning surface             [cross]
benchmarks/loss-decomp/   where nail loses hmmer's hits, by stage [cross]
benchmarks/cutoffs/       per-family score cutoffs from decoys    [reversed]
benchmarks/long-seqs/     how the tools scale with length         [pairs]
benchmarks/pid/           recall against percent identity          [profmark]
benchmarks/thread-scaling/ wall clock against thread count      [cross]
benchmarks/stage-times/   where nail's wall clock goes, by stage  [union]
```

One binary per benchmark, and the shape in brackets is the set it reads.
`cloud-search` and `loss-decomp` name no shape and would read any set with a
query and a target; `cross` is what their labels build. One library sits under
them all, and one binary that belongs to no benchmark: `build-set` makes the
sets.

The warning under **What is tracked** is the one thing in this file to read
before running anything near `sets/`.

### What may be shared

Benchmarks do not share code unless that code is extremely generic. A change
to one benchmark's table, analysis or reader must never be able to break
another. Two benchmarks each carrying their own copy of a table writer is the
right outcome. Do not consolidate it later.

The test is consumers per public item, counted by how many binaries reference
it. An item with one consumer belongs in that consumer. A shared item must
also name no concept belonging to a benchmark: a padded-table renderer, a path
resolver and a manifest reader qualify; anything naming a `pass` string, a
checkpoint or the sensitivities a sweep runs does not.

A tool's parameters belong to the benchmark that sweeps them. Only the
settings that have to be equal for the tools to compare at all stay shared:
one e-value, one thread count, one seeding mode.

`util` is what is left after applying that. Each of `recall`, `cloud-search`
and `loss-decomp` holds a private `src/scores/` module with only what its own
benchmark reaches.

What stays shared is the envelope: where a run's directories are, how a query
is split and its parts joined back up, which binary to run, how a command is
tagged so the ledger can read it back. A tool's own flags belong to the
benchmark that runs it, except where more than one runs it the same way:
`nail`, `hmmsearch` and `mmseqs` have builders in `util::search` for that
reason.

No benchmark looks on `PATH`. `util::tools` holds the path to every binary and
every download, and a benchmark reads it rather than guessing, so a run uses
whichever `hmmsearch` was built for this repo.

### What a benchmark is

A benchmark answers one question. `pid` measures recall against percent
identity; a subcommand that correlates two runs' scores answers a different
question and carries its own reader and output format with it. Work that
answers a different question is a different benchmark, or it is a script, or
it is nothing.

### Structure

The rules the code is being simplified against. We wrote them while taking
`pid` from 2,002 lines to 1,626 and deleting two crates. Auditing the rest of
the repository against them is an open foam issue, and the rules will change
as that audit finds cases they do not settle.

- **Most nouns do not need a type.** Start with a function and its arguments.
  A struct built at one call site and consumed at the next is arguments. A
  struct holding the two things a function returns is a tuple. An enum parsed
  from a string and printed back as the same string is the string. A module
  whose job is to join paths is those joins, written where they run. A struct
  wrapping N calls to a function is those N calls.
- **A doc comment arguing for a type is evidence the type should not exist.**
  `/// passed together because they always travel together`, over three
  `&str` fields, describes three arguments.
- **Abstract what repeats, not what is easy to wrap.** Two search helper
  crates wrapped directory naming, query splitting and binary lookup, and left
  the part written out seven times per benchmark at every call site: tag the
  command, wrap it in a step, push it.
- **Delete dead code rather than annotating it.** An `#[allow(dead_code)]`
  over an item, or a comment saying nothing reads a field, marks something to
  delete.
- **Write a rule once.** Two copies of one calculation drift. One of ours
  carried a `TODO` the other did not.
- **One way per file to do one thing.** A file that opens a result table
  three ways duplicates nothing and is still three shapes to learn.
- **A helper that grows an argument per caller has the disease it was meant
  to cure.** A type per noun and a parameter per caller are the same mistake.
  When `hmmer()` reached ten arguments it was close to being two functions.

### Method

- **Let the compiler find what is dead.** Moving a crate into a binary as a
  module turns its unused public surface into warnings, and those warnings
  are the worklist.
- **Prove a refactor by diffing bytes.** Build a toy set, run it, keep the
  analysis tables. Refactor. Re-run what changed and compare the tables byte
  for byte. Wall-clock columns move; nothing else may.
- **Work in a sandbox copy**, as **Testing behaviour** describes. A refactor
  that reruns a pipeline writes into the store.

### Which build of a tool ran

`make nail` installs the pinned tarball. `make nail NAIL_SRC=/path/to/nail`
builds from a working tree instead, which is how an unreleased fix gets used.
One nail is installed at a time, and `tools/installed.tbl` records what each
rule installed: the version, a sha256 of the binary and where it came from.

The hash is the part that matters. A version string cannot tell two builds
apart, and a tree's binary changes whenever the tree does. So `make check`
compares the hash on disk against the one recorded, and
`util::tools::identity` computes the same digest at run time rather than
reading the table, which is the only way it cannot be stale.

A pipeline stamps what it ran with: `ledger.tbl` gets a `#= tool <name>
<version> <hash>` line per tool when the run finishes, and `parse` carries
those into `runs.tbl` and `scores.tbl`. The run is the only moment this can be
read honestly, since an analysis may happen after a rebuild. A table written
before this existed has no such line and still reads.

## External crates

- `michi` runs a pipeline of commands, times each one, and writes what it did
  to `manifest.tbl`. Its leases are node-aware: a command asking for cores
  gets them off one memory node, and the table says which node and which
  memory policy each command ran under, in `node` and `policy` columns that
  `util::manifest` counts as michi's accounting rather than a run's settings.

  Every michi-batched hmmer here runs in a shared pool: `util::search::hmmer`
  gives its batch `Step::pool(HMMER_CPU x parts)` and its parts no cores of
  their own, so each part runs across the whole pool. A `--cpu 2` hmmsearch
  keeps more than two cpus busy, and on private 2-core leases a batch runs
  1.5x slower. nail and mmseqs are not batched this way; time them before
  switching.

  No core count is written into the code. Every pipeline, build-set's
  included, runs in a michi pool of `--threads` cores, and that value has no
  default: a command that asks for no cores of its own shares the pool, and
  one that does leases out of it. A `default_value_t` on `--threads`, a thread
  count in a `paths.toml` recipe and a fallback to `available_parallelism`
  all count as hard-coding. thread-scaling's pool is its top rung, and
  `--rungs` has no default either. cutoffs' fanout sizes its pool by `--jobs`,
  since its per-family searches run one thread each. The tools' own thread
  counts are another matter: `HMMER_CPU` is 2, and cutoffs' fanout runs each
  search at 1, both on purpose, and say so when adding one. The plot commands
  are the only pipelines left unpinned.
- `libsail` reads and writes the formats: FASTA, Stockholm, p7hmm, and the hit
  tables nail, HMMER, MMseqs2 and BLAST produce. It also draws the samples
  `build-set` deals from: `sample_in_order(m, seed)` yields a uniform draw in
  the order the file holds it, in constant memory. The draw arriving in file
  order is why the deal shuffles which shard each record lands in; a fixed
  round robin over an ascending draw is a stride rather than a partition.
- `toil` writes the padded, `#`-headed tables, and every one of them here is
  written and read through it, michi's manifest included. A `#=` line is a key
  and its words, written with `meta(key, words)` and read back as a `MetaRow`,
  so nothing here splits one by hand. The score tables, too big to hold, are
  read a row at a time through toil's `Reader` in each benchmark's
  `scores/frame.rs`.

## The shape a benchmark has

Three stages, joined by data rather than by code, and one tree holding what
each of them produced. Where an artifact lives follows from what it is, not
from which crate made it, so a set built by one benchmark is readable by
another without either naming the other's directory.

A set and everything derived from it share a directory, so what a run was
searched against is the directory it sits in. The figures are the exception,
and the reason is below:

```
sets/<benchmark>/
├── inputs/
│   ├── set.tbl                  one row per search unit, and what it is made of
│   └── ...                      the queries and targets it names
├── outputs/
│   ├── manifest.tbl             every command, its wall clock, its exit code
│   ├── ledger.tbl               one row per run per shard, and what it cost
│   └── results/<run>.<shard>.tbl   one hit table per run, per target shard
├── analysis/                    what parse worked out
└── tmp/                         scratch, and nothing worth keeping

figures/                         the pdfs, outside the set
```

A benchmark that runs in one pass writes straight into `outputs/`. cutoffs
names a directory under it per stage, `recruit/`, `gather/`, `reject/`,
because each of those is a pipeline of its own with its own manifest and its
own results; pid names `search/` and `reject/`. So the level, where it exists,
is a stage rather than a repeat of the benchmark's name.

The figures are the one thing that does not live with the set. A pdf is read
by a person rather than by another pipeline, so digging five levels down to
find one is cost with nothing on the other side of it. They go to `figures/`
at the root, flat, named for what they draw.

A toy is a whole tree of its own rather than a suffix: `sets-toy/recall/`
beside `sets/recall/`, and `figures-toy/` beside `figures/`. That is what lets
the set directories and the figure names carry no `-toy` or `-real`, and what
keeps a toy figure from overwriting the real one it is named the same as.

Every toy has to work end to end: build, run, parse, plot, with real hits at
every stage. A toy that runs clean and finds nothing leaves everything after
the search untested. Scale the toy, never add a fixture beside it.

Nothing in the code knows that layout. It is what the `paths.toml` files happen
to say, and moving a set to a scratch disk is editing a line rather than
changing anything that compiles.

The scratch sits beside the record rather than inside it, so `outputs/` holds
the record and only the record, and the whole of `tmp/` can go at any time
without touching it. The search pipelines take `--tmp` to put their own
scratch somewhere else, a scratch disk being the usual reason; `build` and
`cutoffs` have no such flag.

Nothing here removes a run. `build-set --rebuild` takes a set back, printing
how many files and how many bytes are about to go and waiting for a `y`, and
what a run left is removed by hand. Never point any of this at recall's `real`
label: see **What is tracked**.

### set.tbl, and why it looks like ledger.tbl

`set.tbl` is what makes a pipeline independent of the recipe that built what it
searches. One row per search unit, one query-and-target pair a tool will be
asked to run, with a fixed spine of `unit`, `query_hmm`, `query_sto`,
`query_fa`, `query_db` and `target`, and every other column an attribute the
recipe wrote down: a shard number, a rung, a family, a residue count. Paths are
relative to the set's own directory, so a set is one tree that can be moved or
linked without rewriting the table. An empty cell means the builder did not
produce that representation, and asking for it fails naming the set rather than
failing later on a path that was never there.

That is the same bargain `ledger.tbl` strikes one seam later, and deliberately
so. A ledger row is a spine plus an open map of settings, which is what lets an
analysis read a run's shape out of the table rather than out of the filenames.
A set row is a spine plus an open map of attributes, which is what lets a search
read a set's shape the same way. `fixed` and `cross` produce the same table
with different attribute columns, and a pipeline reading it does not learn
which one ran.

Both are written by the producer in the pass that produces the artifact, and
nothing else ever writes them. That discipline is the whole defence against a
manifest that says one thing while the directory says another.

`build-set` counts each shard as it deals it and writes the counts into
`set.tbl`: counting a thousand shards afterwards is the whole deal read again,
and the count is only a metadata line.

Before giving anything a `set.tbl`, ask whether `build-set` could produce it
from a label. A set is what a recipe makes from sources, deterministically;
what a search happened to score is never a set.

### The shapes

An open manifest lets one format describe a deal of shards and a grid of
sources. The cost is that a set and the benchmark reading it can disagree with
nothing saying so, and a `reversed` set searched as if it were `fixed` reports
every score against a sequence written backwards. A shape closes that: the
builder stamps `#= shape` and the benchmark names the one it reads, and
`Set::load_as` holds them to each other before a tool runs.

| shape | representations | attributes | read by |
|---|---|---|---|
| `fixed` | `query_hmm`, `query_sto`, `query_db`, `target` | `shard`, `seqs`, `residues`, `bytes` | recall |
| `reversed` | the same as `fixed` | the same as `fixed` | cutoffs |
| `cross` | the same as `fixed` | `query_src`, `target_src`, `seqs`, `residues`, `bytes` | cloud-search, loss-decomp |
| `pairs` | `query_fa`, `target` | `pair`, `query_residues`, `residues` | long-seqs |
| `profmark` | `query_hmm`, `query_sto`, `query_fa`, `target` | `truth`, `originals` | pid |
| `union` | `query_hmm`, `target` | `part` | stage-times |

`reversed` is `fixed` written backwards: the same sources, the same seed and
the same sharding, with each sequence reversed as it is dealt. Reversing keeps
a sequence's composition and destroys its homology, so the two sets hold the
same draw and answer different questions, which is why the columns are
identical and only the declared shape tells them apart. A `fixed` recipe earns
it with `reversed = true`.

`profmark` is one unit: every query against one target file, so the set has a
single row. What tells a true pair from a decoy, and at what percent identity,
is per pair rather than per unit, so it cannot be a column. `truth` names the
file that carries it, relative to the set root, the way the representation
columns name theirs. `originals` names the decoys unreversed, the same way.
The shape check sees that each file is named, not what is in it.

`union` is other sets as one. Each part is a whole recipe of its own, built
under its own directory with its own `set.tbl`, and the union's rows are the
parts' units under the part's name with only the query profile and the target,
since a part's other columns mean what its own shape says. It is for a
benchmark whose units are of different shapes.

`cross` is the only shape where both sides are lists of independent sources
rather than cuts of one. `fixed` grids a query over a partition of a single
target source, so a unit there is a piece of the same thing; in a `cross` a
unit is one whole source against another, and what moves between units is
which sources they are. `pairs` is the diagonal of a cross that was never
built. N x M, N x 1 and 1 x M are one recipe and one manifest, a row per
combination, `unit` naming both sides.

cloud-search and loss-decomp read `cross` because their questions are about
the kind of sequence being searched rather than about how much of it there is.
cloud-search draws Pfam against MGnify and against Swissprot, 500,000 a side,
so its two surfaces differ only in where the targets came from; Swissprot is
not much bigger than that, and a target the size of most of its source is a
sample of nothing. loss-decomp draws Pfam against one recall shard's worth of
MGnify, 2,455,940 sequences, and nothing else.

`pairs` is the one shape with no draw in it. Its two sources are separate
directories of fasta, so a sequence is never on both sides, and pair `i` is the
`i`th file of each in name order. These are pairs somebody chose for their
length, so the recipe selects and places them rather than sampling.

They live in `util::set::shape` rather than with a benchmark, because a shape is
the agreement between a builder and a reader and neither side owns it. `fixed`
carries `seqs`/`residues`/`bytes` because the analyses read them: a shape is
what a whole benchmark needs, run and parse together, not what its first stage
opens.

A set built before shapes existed says nothing about itself, and is then checked
on its columns alone.

Two benchmarks name no shape at all. cloud-search and loss-decomp read
`query_hmm` and `target` and no attribute, so `Set::load_needing` holds them to
those and to nothing else, and either would read any set that has them. Their
own labels build a `cross`. A reader that wants an attribute names the shape
instead, because then the recipe is what it depends on.

## paths.toml

No crate resolves a location. A tool is told where its inputs are and where its
outputs go, and the telling is a `paths.toml` in the tool's own crate directory,
one table per label:

```toml
[toy]
set      = "../../sets-toy/cloud-search/inputs"
run      = "../../sets-toy/cloud-search/outputs"
analysis = "../../sets-toy/cloud-search/analysis"
figures  = "../../figures-toy"
tmp      = "../../sets-toy/cloud-search/tmp"
```

A benchmark that draws nothing has no `figures` key and no field for it:
recall, cutoffs and long-seqs stop at `analysis`.

A label is a whole set of paths under one name, so a toy run and a real run
differ by a word: `recall run --in toy --threads 8`. Running a tool without
`--in` prints the labels its file holds.

Every benchmark has `toy` and `real`, over two sets of its own:
`sets-toy/<benchmark>/` and `sets/<benchmark>/`, which `build-set` builds
under the labels `<benchmark>-toy` and `<benchmark>-real`. cutoffs has four
more, `sp-10k` through `sp-all`, which are neither toy nor real and sit in
`sets/` beside the rest.

No benchmark reads a set another one built.

Relative paths resolve against the file's own directory, which is why a crate
needs no notion of a repository, and an absolute path is left alone: a set on
a scratch disk is one line.

**What may go in one of these: paths, and for a tool that makes a dataset, how
much of it to make.** `build-set`'s labels carry `shards`, `seqs` and rungs,
because a toy set is defined by being small and a size is not a path. Nothing
about a search belongs in any of them: no tools, no flags, no sensitivities,
no threads, no name templates. A strategy such as cutoffs' `--strategy` is
CLI only.

That line is drawn where it is because of what the last generation of these
files did. They built argv out of templates, which made the file the thing
that decided what ran, turned the Rust into an interpreter for it, and left
the real logic in strings nothing type-checked. Every value read now is
deserialized into a typed field with `deny_unknown_fields`, so a key that does
not belong is a parse error naming the label rather than something silently
ignored.

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
since nothing here ran those commands and an exit code or an argv would have
to be invented. They arrive with `.time` files beside the tables, and the
deleted `store import` read those into a ledger; it was the only thing that
wrote a ledger by hand.

## benchmarks/util

- `cut` cuts an hmm file and its alignments by name: `subset_*` into one file,
  `scatter_*` into one file per record. The sibling of `split`, which cuts the
  same files by weight.
- `paths` reads the `paths.toml` a tool keeps beside its own source: the labels
  it can run under, and what each names as input and output. Nothing else in
  this workspace resolves a location.
- `set` holds the shape a pipeline reads: one row per search unit, with the
  query, the target and whatever the recipe wrote down about each. It reads and
  writes `set.tbl`.
- `manifest` reads back the table `michi`'s sink wrote, and builds the
  `results/` paths from it.
- `ledger` holds the shape the analyses read: one row per run per shard, with
  each run's seconds, core-seconds and peak resident set already totalled. It
  distills a manifest into that shape, and reads and writes `ledger.tbl`. The
  three fold differently. Batched commands overlap, so the step takes as long
  as its slowest and holds every one of their resident sets at once, while
  serial commands take their total and only ever hold one. Core-seconds are
  work rather than elapsed time, so they add however the commands were
  scheduled.
- `tools` holds where the binaries and the downloads are.
- `split` cuts a query set into balanced parts for a batch of jobs.
- `search` builds the commands a pipeline is assembled out of: the directories
  a run writes into, the split and the cat that put a query's parts back
  together, and a builder each for nail, hmmsearch and mmseqs. It also holds
  the four settings every benchmark is held to so that its tools compare,
  `SEED_MODE`, `SEED_S`, `HMMER_CPU` and `EVALUE`, and nothing else about a
  sweep. A tool's parameters belong to the benchmark that sweeps them.
- `clean` measures a list of directories, shows what is in them, asks, and
  deletes.
- `domains` classifies how a profile's HMMER domain hits on one sequence make
  up its score: one domain, overlapping copies of one model region, domains in
  order, or domains out of order. `measure` takes no cutoff and keeps every
  candidate score; `classify` holds a measure against a per-family cutoff and
  returns the May 2026 analysis's seven categories exactly; `pattern` names
  the arrangement that scores most, for results with no cutoff, and adds
  `diffuse` for a pair with no domain of 4 bits or more.

## The Pfam-against-MGnify benchmarks

Pfam profiles against MGnify metagenomic sequences, and the largest sets here
by a long way. `build-set` cuts the sources into a set under `sets/`: `fixed`
is one query set against target shards of equal size, `cross` is every query
source against every target source. Each writes a `set.tbl`, and every
pipeline here takes `--in <label>` to say which one to search.

Three benchmarks. `recall` searches every shard while sweeping nail's and
MMseqs2's prefilter sensitivity. `cloud-search` seeds once, then searches every
`(A, B)` pruning cell off those same seeds, so the pruning parameters are the
only thing moving. The grid is run twice, at `-a 5` and at `-a 1`, because a
surface at one `-a` cannot say whether a hard-pruning cell found its hits or
recovered them. `loss-decomp` asks where nail loses the hits HMMER finds,
and how far down the prefilter list the hits it does find sit. Its arms are
sensitivities, `--s 12.0,10.0,7.5`, each seeded static with
`--mmseqs-max-seqs` unbounded, so an arm aligns everything the prefilter
returned and its seed list is the most nail could get at that sensitivity.
An arm is one `nail search`, seeding and alignment together, with the seed
list written beside its table so `parse` can tell a pair the seeding never
offered from one it offered and nail dropped; hmmer runs once outside the
sweep, because the truth set is the same for all of them and is most of the
wall clock. An arm leaves mmseqs' databases under
`results/prefilter.<arm>.<unit>/`, which is part of the record: `parse depth`
reads a pair's rank in its query's prefilter list out of it.

There are two table grammars, and which one a pipeline gets is settled by
whether it moves one tool's parameters in a way that changes what that tool
scores a pair. Each of the three benchmarks carries its own copy of the code
that writes and reads them, under its own `src/scores/`, holding only the half
it uses: recall has `write.rs` and `read.rs`, the two sweeps have `runs.rs`.

`recall parse scores --in <label>` reads recall into `scores.tbl`: one row per
query/target pair, one score column per tool, and a `pass` string holding one
character per run. A prefilter sweep changes which pairs a tool reports, not
what it scores them, so one column per tool is the honest shape and the cheap
one; a column per run is what made an earlier grammar unwritable at a
thousand shards.

`cloud-search parse runs --in <label>` reads cloud-search and loss-decomp into
`runs.tbl`: one score column per run, and no pass string. A cell is the score,
or `-` where that run's seeding offered the pair and the run did not report
it, or `.` where the seeding never offered it at all. Whether a run cleared its
family's cutoff is its score against the cutoff `#= cutoffs` names, worked out
when the table is read. A seeding sweep gives every arm its own seed list, so
the offered-or-not answer differs between runs of one pair; a run's ledger row
names the seeding it replayed in a `seeds` setting, the list is
`results/seeds.<seeding>.<shard>`, and many runs to one seeding is the
ordinary case. `-A` and `-B` constrain the dynamic programming, so two cells
can score one pair differently, and that is what the sweep measures. These
pipelines search one shard, so the wider row costs nothing.

A run's name is stacked over its column on the dashes it is built from, so
`A10.0-B12.0-a1` takes three header lines rather than fourteen columns of
padding on every row.

Both are collected the same way: each shard on its own, written as a block in
shard order, so the file comes out the same bytes however many threads it ran.
`--threads` and `--mem` size that. `summary` asks only whether each run
cleared its family's cutoff and what the domain list says, which recall's copy
reads out of the pass string and the sweeps' out of the score columns.
`stages` reads `runs.tbl` alone, because where a pair was dropped is a
question about a run rather than about a tool, and it counts a pair as found
only if nail kept it over the family's cutoff. `summary` still counts what a
run reported; whether it should follow is an open issue.

    recall        parse scores → parse summary
    cloud-search  parse runs   → parse summary → plot
                               → parse stages
    loss-decomp   parse runs   → parse stages  → plot
                               → parse depth

`plot` draws one figure off cloud-search's summary: the `(A, B)` heatmaps of
sensitivity and wall clock, with the `full` run in the far corner. `--full-dp`
records no `-A` and no `-B`, so it has no cell; the corner is where the surface
is heading as the pruning relaxes, and it is the ceiling the pruned cells are
read against.

Each field is drawn at each `-a` and then as the difference between them, six
panels. The difference is what the second arm is run for: it says where on the
grid the sensitivity credited to `-A` and `-B` was recovered by retrying
disjoint clouds instead. A summary with one `-a` in it draws two panels.

`depth` bins each pair's rank in its query's prefilter list in ranges that
double from `--bin`, the way prog's `n_take` does, and writes `depth.tbl`: per
unit, run and bin, the prefilter pairs at that depth, the seeds, the pairs
nail kept over the cutoff, the pairs hmmer found, and the hits, kept and found
both. Hits over prefilter pairs is what a stopping rule is betting on. A last
row per run counts hmmer's pairs the prefilter never returned.

`stages` writes one row per (unit, run) with a column per checkpoint: not in
the mmseqs prefilter list, in it but dropped by mmseqs' alignment, seeded
but dropped by nail's cloud search or forward filter, and scored under the
family's cutoff. The first needs each arm's prefilter database, which is why
`parse stages` takes the run directory as well as the table. It writes the
table three times: `stages.tbl` over every hmmer hit, `stages-single.tbl`
over the hits hmmer reports exactly one domain for, and `stages-multi.tbl`
over the rest. The figure reads the first alone. That is the shape a figure reads. loss-decomp's figure shows where in nail's
pipeline the hits were lost; a point per arm that gives the total without
saying which stage took it is not that figure.

## benchmarks/cutoffs

`cutoffs` calibrates scores. What makes it a calibration rather than a
benchmark is that its product, `data/mgy-cutoffs.tbl`, is committed and
promoted by hand.

### The words

Naming here is deliberate and worth holding to, because the loose version of it
hides the one thing the method turns on.

- **initial reversals**: the source sequences written backwards. These are
  not decoys. Nothing has been asked about them yet.
- **recruits**: initial reversals that scored against some family in the
  first search. `candidate` is not the word: this repo already uses that for
  what a prefilter promotes.
- **decoys**: recruits whose *original* has no significant match to the
  family that recruited them. Unlikely to be homologs, so they are noise, and
  noise is what a threshold is learned from.
- **rejects**: recruits whose original *does* match. Reversed homologs
  wearing a decoy's clothes, and they are thrown out.
- **original**: a recruit re-reversed, which is bit-identical to the source
  sequence. Never *forward*: nail's Forward algorithm owns that word here, and
  a reversal and its original are two forms of one sequence rather than two
  directions of anything. The manifest column is `form`.

### The stages

`recruit` searches every family against the initial reversals, cheaply, and
finds the small subset that scores at all. `gather` pulls those out of the
shards in both forms. `reject` searches each family against its recruits with
the prefilter effectively off, so a score is a real score, and the two forms
together split the recruits into decoys and rejects. `learn` turns the decoy
scores into the per-family cutoffs every hit is afterwards held against.

The set arrives reversed, built by a `fixed` recipe under a `reversed` tag, so
no stage here reverses anything and no second copy of the shards is made. A
record read out of a reversed shard is already the reversal, and reversing it
again gives the original, so both forms come out of one pass in `gather`.
Everything the calibration makes is an output of that set.

The recruits are not a set and get no `set.tbl`: `build-set` makes a set from
sources under a recipe, deterministically, and what `recruit` happened to
score could never be a recipe.

### Rejection is a per-pair test, and that is not negotiable

**A family's null holds only sequences vetted for that family. Never pool them
across families.**

`reject` asks exactly one question about one pair: family A recruited sequence
S, so does A also match S's original? A yes makes S a reject, the reversal of
a true member of A, rather than a decoy for A. Family B's answer for B's own
recruit settles it for B and for nothing else.

The reason this matters more than it looks: a reversed sequence keeps a
surprisingly high similarity score against its unreversed self. That is not
only biology; a good part of it is approximate-palindrome statistics of text,
which the Wheeler lab wrote up in *wasitamatchisaw*. So the reversals of a
family's own true members are exactly the sequences that score well against it
on the reverse pass. The decoys that look best are the ones that are not decoys
at all.

Searching all of Pfam against the union of every family's recruits in one
invocation is fine, and it is how the fanout is avoided. What is not fine is
letting the union widen any family's null. `learn` joins each hit back to the
family that recruited the sequence, and keeps only those.

**The restriction costs nothing, and it comes with its own test.** Widening the
pool cannot change a family's null: the null is the scores of sequences the
family *hit*, and a sequence A never hits contributes nothing, whatever else
sits in the target file. So pooling should move no cutoff at all.

If pooling ever does move a cutoff, that is not an improvement and not noise.
It means pairs reached A's null without being asked A's question, and by the
paragraph above the pairs most likely to do that are the reversals of A's own
true members. The null would then be seeded with disguised true positives, A's
cutoff would rise, and real hits would be discarded quietly, the exact failure
rejection exists to prevent. A measurable effect from pooling is the method
reporting its own unsoundness, so treat it as a bug to find rather than a
result to keep.

The one case where B's recruits could legitimately inform A is A and B being
highly related, and that is the two families not being independent rather than
an argument for pooling.

### `--strategy fanout|union`

Two ways for `reject` to score the recruits. `fanout` is one search per
family per form. `union` is two invocations, all of Pfam against every
recruit's reversal and then against every recruit's original, with `learn`
joining each hit back to the family that recruited the sequence. The per-pair
rule above is what `union` must not break.

Which is faster is unmeasured on MGnify, and nail's and mmseqs' runtimes do
not scale linearly in target size, so do not pick between the arms from an
extrapolation or from the SwissProt result: measure both on the same input.
`auto` needs a fitted model before it means anything, and stays unimplemented
until then.

## benchmarks/pid

Recall as a function of the percent identity between a query and its target.
`build-set --in pid-toy|pid-real` assembles it from a profmark split: Pfam
families divided by identity, with their true targets hidden among whole
TrEMBL sequences written backwards, 100 per true pair. The decoys are drawn so
that their lengths follow the true targets' histogram in 50-residue bins:
drawn uniformly, TrEMBL runs far longer than a Pfam domain, and every tool but
nail's sequence mode ranked the longest decoys highest. `truth.tbl` records
which pair is which and at what identity, and `originals.fa` holds each decoy
unreversed under the name `target.fa` gives it.

The whole benchmark is one search unit, every query against one target file,
so the set has one row, and the truth is per pair rather than per unit. That is
why `truth` names a file instead of being a column: see the `profmark` shape.

The split itself sits at `benchmarks/pid/profmark/` rather than in the set. It
depends only on the alignments and the split parameters, both labels draw from
the same one, and the recipe points at it. It is not in a fresh checkout;
`build-set` draws one when it is missing.

`pid run` searches every tool against the set, into `outputs/search/`. `pid
reject` then searches the family profiles against the originals of the decoys
the runs have to settle, into `outputs/reject/`. `parse` turns the results
into the tables the plot scripts read, and `plot` draws them. All four take
`--in <label>`. pid builds its own commands for the four tools nobody else
runs, blastp, psiblast, lastal and diamond, and takes nail, hmmsearch and
mmseqs from `util::search`.

Rejection is the cutoffs rule applied per pair, with one judge for every
tool: a (query, decoy) hit is dropped when the query's family profile hits
that decoy's original at E ≤ 1e-3 in hmmsearch. A sequence query `family|id`
is judged by its family's profile. The originals searched are the union, over
every run, of the decoys ranked at or above that run's worst true pair. A
decoy ranked below every true pair cannot move a point of the ROC, so the
curve means what it did before rejection existed. The judge runs at `-Z`
equal to `target.fa`'s sequence count, so its E-values are on the scale of the
search they settle rather than of the smaller file it reads.

`parse` counts a false positive once for every true pair of the query that
found it, so a profile query's decoy hit fills as many slots as its family has
true pairs, up to 10. This is the pre-refactor design and stays.

`parse` has three subcommands. `recall` writes `pid.tbl`, `roc.tbl` and
`time.tbl`: one row per run per identity bin, one row per run per point of its
ROC curve, and one row per run of what it cost against what it found. `score`
and `cells` answer questions other than recall against identity.

`plot` names five figures, roc, pid, time, cells and score, and `--only`
draws a subset, one name per flag.

## benchmarks/long-seqs

How each tool's runtime scales as sequences get longer. Six paired
query/target files, where `run` searches each query against its pair and
`parse` turns that into tables. Its inputs are small and checked in under
`data/long-seqs/`, and `build-set --in long-seqs-toy` or `--in
long-seqs-real` cuts them into a set the same way every other benchmark gets
one. Its queries are sequences rather than profiles, which is what the
manifest's `query_fa` column is for, and its shape is `pairs`: one query
against one target, with the two sources separate directories so a sequence
is never on both sides.

## benchmarks/thread-scaling

Strong scaling for nail, mmseqs and hmmer over its own `cross` set, all of
Pfam against 500,000 MGnify sequences, built by `build-set --in
thread-scaling-toy|thread-scaling-real`. `run` searches every unit at every
rung of `--rungs`, `--reps` times over, with the reps as the outer loop so
load from other users spreads across rungs. Four arms: nail `-t N`, mmseqs
`--threads N`, one `hmmsearch --cpu N`, and `hmmer-split`, the query cut N/2
ways at `--cpu 2`. Each search is pinned to N cores. `/proc/loadavg` is read
either side of every search into `outputs/load.tbl`, and it counts the
search's own threads as well as anyone else's.

Every arm of a scaling curve starts from the same rung, one thread. Where an
arm's recipe degenerates at a rung into another arm's command, reuse that
measurement rather than leaving the rung out: `run` gives the split arm no
one-thread search, because one part at `--cpu 1` is `hmmer-t1`, and `parse`
takes `hmmer-t1` as its first point.

`parse` writes `scaling.tbl`, one row per unit per arm per rung with median
wall and core-seconds, peak RSS, the time perfect scaling from the arm's
lowest rung would take and the percentage of that speedup reached, and the
highest load seen; and `agree.tbl`, one row per run counting the pairs
missing, extra or rescored against rep 1 of the arm's lowest rung. `plot`
draws wall clock and that percentage against threads into
`thread-scaling.pdf`.

The split arm is `util::search::hmmer`, so its parts share one pool of N
cores.

## benchmarks/stage-times

Where nail's wall clock goes, stage by stage, and how that changes with the
input. One nail search per unit at the shared settings, with `-s` on, and no
sweep: the units are the inputs. Its set is a `union` of two parts: `mgy`, a
`cross` of all of Pfam against one recall shard's worth of MGnify, 2,455,940
sequences, and `pid`, the profmark set pid searches, drawn from the same split
with the same parameters. Both run at 48 threads.

`-s` makes nail print a stage tree to stdout when it finishes: `setup`,
`seeding` and `alignment`, each in wall clock, and under each the stages as
leaves. Under `align` in prog mode it prints one line per round. michi
discards stdout unless told where to put it, so `run` sends it to a `.stats`
file beside the unit's hit table.

The seeding and setup leaves are wall, since each is one serial step or one
mmseqs invocation. The alignment leaves run on every thread, so nail prints
each one's cpu seconds, the sum over every pair across every thread, beside a
wall that is its share of that sum times the branch's wall. That share-times-
wall is exact when every thread is busy for the whole branch, and the branch
line says how busy they were. `parse` writes `stages.tbl`, one row per unit
per line of the tree: its parent, its wall, its share of its parent, and its
cpu where nail printed one. `counts.tbl` carries the counts nail prints above
the tree, and a `#= michi` line per unit carries what michi timed for the
whole command. A `.stats` file from a nail before the three-branch tree reads
wrong, since its alignment leaves were cpu seconds; rerun rather than reparse.

## Formatting

rustfmt is the standard here. Run `cargo fmt --all` and commit what it does.

A file it reformats is a file that was committed unformatted, so the churn is
a fix rather than noise, and keeping it out of a commit to make that commit
read better only leaves the next person to do it. There is no house style that
overrides it and no file exempt from it.

## Testing behaviour: work in your own copy

This working tree belongs to whoever is at the keyboard. A pipeline writes into
`sets/`, a `build` subcommand writes a set there, and when two people write
there at once neither can tell which results are theirs. So run nothing here.
Copy the source into `tmp-claude/sandbox/`, link the expensive directories, and
work in the copy.

```bash
ROOT=$(git rev-parse --show-toplevel)
SB=$ROOT/tmp-claude/sandbox

rsync -a --delete \
  --exclude '.git/' --exclude 'target/' --exclude 'tmp-claude/' \
  --exclude '/data' --exclude '/tools' --exclude '/benchmarks/pid/profmark' \
  --exclude '/sets' --exclude '/sets-toy' \
  --exclude 'outputs/' --exclude 'tmp/' \
  "$ROOT/" "$SB/"

ln -sfn "$ROOT/data" "$SB/data"
ln -sfn "$ROOT/tools" "$SB/tools"
# the split is drawn lazily, so the real tree may not have one, and a
# sandbox that drew its own would get the link nested inside it
if [ -d "$ROOT/benchmarks/pid/profmark" ] &&
   { [ -L "$SB/benchmarks/pid/profmark" ] || [ ! -e "$SB/benchmarks/pid/profmark" ]; }; then
  ln -sfn "$ROOT/benchmarks/pid/profmark" "$SB/benchmarks/pid/profmark"
fi

```

Run that before testing anything, and again after every edit: a copy goes
stale, and a result from stale source is worth nothing. `--delete` is what
keeps it current, and it drops what was deleted from the source without
touching the sandbox's own `sets/`, `target/`, pid `inputs/` or links, since
rsync leaves excluded paths on the receiving side alone. The link paths and
`/sets` are excluded without a trailing slash on purpose: a pattern ending in
`/` matches only a directory, and on the receiving side three of these are
symlinks, so `--delete` removes them. It copies uncommitted edits, which is the
point: what wants testing is usually not committed yet.

The sandbox gets its own `sets/`, which is what makes a build or a run in the
copy safe: `util::tools::repo()` resolves from `util`'s own
`CARGO_MANIFEST_DIR`, so a build inside the copy resolves the copy's root and
writes the copy's sets.

Then run inside `$SB`, through its own shims. Every binary links to the same
`benchmarks/shim`, which builds the crate the link is named after and runs what
it built: `benchmarks/recall/recall run --in toy --threads 8`, and so on. Every
path these crates resolve comes from a `CARGO_MANIFEST_DIR`, so a build in the
sandbox reads the sandbox's `data/` and `tools/` links and writes the
sandbox's `sets/`.

The links are the two big directories plus the profmark split, which is
expensive for the reasons the pid section gives. The sandbox reads all three;
nothing in it should write them, so do not run `make data` or `make tools`
from the copy.

`/sets` and `/sets-toy` are excluded because a build writes them and nothing
commits them, so the copy builds its own. Without the exclude, `--delete`
removes the set a `build-set` in the sandbox just wrote, since the real tree
has nothing there to match it. long-seqs' set is checked in under
`data/long-seqs/`, which is a link, so the copy reads the real one.

`/sets` is the exclude that matters most, and dropping it is expensive rather
than wrong: a checkout that has recall's results holds terabytes there, and
rsync will copy every byte of it into the sandbox. Check this list against
what is on disk before running the recipe.

Edit in the real tree and re-sync, never in the sandbox: an edit in the copy is
gone at the next sync. To check an analysis against a run that finished
elsewhere, copy that run's `sets/<benchmark>/outputs/` into the sandbox and
parse it there.

Prove a refactor by diffing bytes: build a toy set, run it, keep the analysis
tables, refactor, re-run what changed and compare the tables byte for byte.
Wall-clock columns move; nothing else may.

## What is tracked

The downloads under `data/` and the builds under `tools/bin/` are not, and
neither is anything under `sets/` or `sets-toy/`: a set, a run and an analysis
are usually things this repository can make again.

**`sets/recall/` is the exception, and it is not backed up by being
reproducible.** It holds about 3.5 TB of results searched on a cluster, the
hmmer set alone some 21,926 hours of wall clock, and nothing here can produce
them again. Git ignores it like every other set, so the only thing protecting
it is that nobody deletes it. Read it, never write it, and never point an
`rsync --delete` or a `build-set --rebuild` at that label. The deal predates
`set.tbl`, so the set cannot be rebuilt from its own manifest either.

Two exceptions are committed on purpose: `data/long-seqs/`, because nothing
fetches it and ignoring it would empty a fresh clone, and
`data/mgy-cutoffs.tbl`, because learning it is a calibration run rather than a
download.

Plotting is Python and matplotlib. Four benchmarks have a `plot` command and
a `scripts/` beside their source: cloud-search, loss-decomp, pid and
thread-scaling. recall, cutoffs and long-seqs stop at `parse`.

`reports/` is untracked too. Each report there was put together by hand from
`sets/` and scratch analyses, and no command here regenerates one.
