# nail-benchmarks

Benchmarks for [nail](https://github.com/travisWheelerLab/nail), a profile HMM
search tool, against the tools it is compared to: HMMER, MMseqs2, BLAST, LAST
and DIAMOND.

A Cargo workspace at the repo root. Binaries are the interface. There is one
per benchmark, and what each of them builds, runs and works out goes under
`sets/`.

Each binary has a shim beside its `Cargo.toml`, named after it, so the
`recall run` and `pid parse` spellings used throughout this file are what you
type:

```
benchmarks/recall/recall run --in real
benchmarks/build-set/build-set --in recall-toy
benchmarks/pid/pid parse recall
```

Every one of them is a symlink to `benchmarks/shim`; the crate it builds is
the name on the link. It builds that crate and runs it only if the build
succeeded, so a broken build never falls through to the last binary that
worked. Nothing can
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
Makefile                  downloads data, builds tools, nothing else
data/                     what the Makefile downloaded
tools/bin/                what the Makefile built
sets/                     what a build, a run or an analysis produced
sets-toy/                 the same, for the toy of each
figures/                  the pdfs, flat
figures-toy/              the same, for the toys
benchmarks/shim           the build-and-run shim every binary links to

benchmarks/util/          set, search, paths, manifest, ledger, tbl, tools,
                          split, cut, clean

benchmarks/build-set/     cuts sources into a set: fixed | reversed | cross |
                          pairs | profmark

benchmarks/recall/        recall against prefilter sensitivity   [fixed]
benchmarks/cloud-search/  the (A, B) pruning surface             [cross]
benchmarks/loss-decomp/   where nail loses hmmer's hits, by stage [cross]
benchmarks/cutoffs/       per-family score cutoffs from decoys    [reversed]
benchmarks/long-seqs/     how the tools scale with length         [pairs]
benchmarks/pid/           recall against percent identity          [profmark]
```

One binary per benchmark, and the shape in brackets is the set it reads.
`cloud-search` and `loss-decomp` name no shape and would read any set with a
query and a target; `cross` is what their labels build. One library sits under
them all, and one binary that belongs to no benchmark: `build-set` makes the
sets.

The warning under **What is tracked** is the one thing in this file to read
before running anything near `sets/`.

### How much of the libraries is actually shared

Counted per public item, by how many binaries reference it. A part with one
consumer belongs in that consumer rather than in a library, and the count is
the test rather than how well the library reads from the inside.

`util` is what is left after applying that. The `scores` crate was three
benchmarks' tables sharing an envelope, and it is now a private `src/scores/`
module inside each of `recall`, `cloud-search` and `loss-decomp`, each holding
only what its own benchmark reaches. The `search` crate became
`util::search`, and the sweep constants that were in it went to the benchmarks
that sweep them.

What stays shared is the envelope: where a run's directories are, how a query
is split and its parts joined back up, which binary to run, how a command is
tagged so the ledger can read it back. A tool's own flags belong to the
benchmark that runs it, except where more than one runs it the same way --
`nail`, `hmmsearch` and `mmseqs` have builders in `util::search` for that
reason.

No benchmark looks on `PATH`. `util::tools` holds the path to every binary and
every download, and a benchmark reads it rather than guessing, so a run uses
whichever `hmmsearch` was built for this repo.

### Which build of a tool ran

`make nail` installs the pinned tarball. `make nail NAIL_SRC=/path/to/nail`
builds from a working tree instead, which is how an unreleased fix gets used.
One nail is installed at a time, and `tools/installed.tbl` records what each
rule installed: the version, a sha256 of the binary and where it came from.

The hash is the part that matters. A version string cannot tell two builds
apart -- a release and a working tree can both say `nail 0.7.1` and be different
binaries -- and a tree's binary changes whenever the tree does. So `make check`
compares the hash on disk against the one recorded and says when they differ,
and `util::tools::identity` computes the same digest at run time rather than
reading the table, which is the only way it cannot be stale.

A pipeline stamps what it ran with: `ledger.tbl` gets a `#= tool <name>
<version> <hash>` line per tool when the run finishes, and `parse` carries
those into `runs.tbl` and `scores.tbl`. The run is the only moment this can be
read honestly, since an analysis may happen after a rebuild. A table written
before this existed has no such line and still reads.

`pid` reads a set like everything else. What kept it out was a record-level
truth table, and the answer was for `set.tbl` to name the file rather than
carry it -- see the `profmark` shape. It builds its own commands for the four
tools nobody else runs -- blastp, psiblast, lastal and diamond -- and takes
nail, hmmsearch and mmseqs from `util::search`. Its `profmark/` split lives at
the crate root and is not in a fresh checkout; `build-set` draws one when it is
missing.

## External crates

- `michi` runs a pipeline of commands, times each one, and writes what it did to
  `manifest.tbl`. It pins each command with `sched_setaffinity` and sets no
  memory policy, handing out one logical cpu per physical core from the low end
  of the pool. On this two-node box that lease straddles both nodes, since node
  membership goes by parity. Measured, four reps of a 350s nail search at 32
  cores: straddling costs about 2% against 32 cores of a single node, and
  `numactl --membind` on top of single-node pinning buys nothing (+0.35%,
  inside the noise). So if this is ever worth fixing, the fix is making the
  lease node-aware rather than adding `set_mempolicy`. It is not worth fixing
  today. Full write-up while it lasts: `tmp-claude/numaprobe/REPORT.md`.
- `libsail` reads and writes the formats: FASTA, Stockholm, p7hmm, and the hit
  tables nail, HMMER, MMseqs2 and BLAST produce. Since 0.4.0 it also draws the
  samples `build-set` deals from: `sample_in_order(m, seed)` yields a uniform
  draw in the order the file holds it, in constant memory, where the 0.3 line
  built a vector of every index first -- 19.6 GB of them for MGnify's 2.4
  billion records. The draw arriving in file order is why the deal shuffles
  which shard each record lands in; a fixed round robin over an ascending draw
  is a stride rather than a partition.
- `tabl` writes the padded, `#`-headed tables. It is not published: the
  workspace takes it as a path dependency on a sibling checkout, `../tabl`.
- `feisty` sits in `[workspace.dependencies]` and no member uses it.

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
names a directory under it per stage -- `recruit/`, `gather/`, `reject/` --
because each of those is a pipeline of its own with its own manifest and its
own results. So the level, where it exists, is a stage rather than a repeat of
the benchmark's name.

The figures are the one thing that does not live with the set. A pdf is read
by a person rather than by another pipeline, so digging five levels down to
find one is cost with nothing on the other side of it. They go to `figures/`
at the root, flat, named for what they draw.

A toy is a whole tree of its own rather than a suffix: `sets-toy/recall/`
beside `sets/recall/`, and `figures-toy/` beside `figures/`. That is what lets
the set directories and the figure names drop the `-toy` and `-real` they used
to carry, and what keeps a toy figure from overwriting the real one it is
named the same as.

Nothing in the code knows that layout. It is what the `paths.toml` files happen
to say, and moving a set to a scratch disk is editing a line rather than
changing anything that compiles.

The scratch sits beside the record rather than inside it, so
`outputs/` holds the record and only the record, and the whole of `tmp/`
can go at any time without touching it. The search pipelines still
take `--tmp` to put their own scratch somewhere else, a scratch disk being the
usual reason; `build` and `cutoffs` have no such flag and never had one.

Nothing here removes a run. `build-set --rebuild` takes a set back, printing
how many files and how many bytes are about to go and waiting for a `y`, and
what a run left is removed by hand. There used to be a `store clean` that did
the second; git has it. Never point any of this at recall's `real` label: see
**What is tracked**.

### set.tbl, and why it looks like ledger.tbl

`set.tbl` is what makes a pipeline independent of the recipe that built what it
searches. One row per search unit -- one query-and-target pair a tool will be
asked to run -- with a fixed spine of `unit`, `query_hmm`, `query_sto`,
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

`reversed` is `fixed` written backwards: the same sources, the same seed and
the same sharding, with each sequence reversed as it is dealt. Reversing keeps
a sequence's composition and destroys its homology, so the two sets hold the
same draw and answer different questions, which is why the columns are
identical and only the declared shape tells them apart. A `fixed` recipe earns
it with `reversed = true`.

`profmark` is one unit: every query against one target file, so the set has a
single row. What tells a true pair from a decoy, and at what percent identity,
is per pair rather than per unit, so it cannot be a column -- `truth` names the
file that carries it, relative to the set root, the way the representation
columns name theirs. `originals` names the decoys unreversed, the same way.
The shape check sees that each file is named, not what is in it.

`cross` is the only shape where both sides are lists of independent sources
rather than cuts of one. `fixed` grids a query over a partition of a single
target source, so a unit there is a piece of the same thing; in a `cross` a
unit is one whole source against another, and what moves between units is
which sources they are.
`pairs` is the diagonal of a cross that was never built. N x M, N x 1 and 1 x M
are one recipe and one manifest, a row per combination, `unit` naming both
sides.

cloud-search and loss-decomp read it because their questions are about the kind
of sequence being searched rather than about how much of it there is: both draw
Pfam against MGnify and against Swissprot, 500,000 a side, so their two surfaces
differ only in where the targets came from. Swissprot holds 570,420 sequences
in total, which is what sets that size -- a target the size of most of its
source is a sample of nothing.

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

`util::set::shape` also declares `ladder`, nested rungs on both axes. Nothing
builds one and nothing reads one since `calibrate` went, so it and its recipe
in `build-set` are dead code.

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
differ by a word: `recall run --in toy`. Running a tool without `--in` prints
the labels its file holds.

Every benchmark has `toy` and `real`, over two sets of its own:
`sets-toy/<benchmark>/` and `sets/<benchmark>/`, which `build-set` builds
under the labels `<benchmark>-toy` and `<benchmark>-real`. cutoffs has four
more, `sp-10k` through `sp-all`, which are neither toy nor real and sit in
`sets/` beside the rest.

No benchmark reads a set another one built. cloud-search used to, running its
grid over the profmark split pid assembles, and that label is gone.

Relative paths resolve against the file's own directory, which is why a crate
needs no notion of a repository, and an absolute path is left alone -- a set on
a scratch disk is one line.

**What may go in one of these: paths, and for a tool that makes a dataset, how
much of it to make.** `build-set`'s labels carry `shards`, `seqs` and rungs,
because a toy set is defined by being small and a size is not a path. Nothing
about a search belongs in any of them -- no tools, no flags, no sensitivities,
no threads, no name templates.

That line is drawn where it is because of what the last generation of these
files did. They built argv out of templates:

```toml
args = "--allow-overwrite --mmseqs-s {s} --seed-mode prog -E {evalue}"
```

which made the file the thing that decided what ran, turned the Rust into an
interpreter for it, and left the real logic in strings nothing type-checked.
Every value read now is deserialized into a typed field with
`deny_unknown_fields`, so a key that does not belong is a parse error naming the
label rather than something silently ignored.

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
to be invented. They arrive with `.time` files beside the tables, and a
`store import` used to read those into a ledger. It was the only thing that
wrote a ledger by hand, and `sets/recall/outputs/ledger.tbl` is what it wrote.
The crate is deleted; restore it from git when a cluster run next lands.

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
  serial commands take their total and only ever hold one. Core-seconds are work rather than
  elapsed time, so they add however the commands were scheduled.
- `tbl` writes the padded, `#`-headed table every analysis produces.
- `tools` holds where the binaries and the downloads are.
- `split` cuts a query set into balanced parts for a batch of jobs.
- `search` builds the commands a pipeline is assembled out of: the directories
  a run writes into, the split and the cat that put a query's parts back
  together, and a builder each for nail, hmmsearch and mmseqs. It also holds
  the four settings every benchmark is held to so that its tools compare --
  `SEED_MODE`, `SEED_S`, `HMMER_CPU` and `EVALUE` -- and nothing else about a
  sweep.
- `clean` measures a list of directories, shows what is in them, asks, and
  deletes.

## The Pfam-against-MGnify benchmarks

Pfam profiles against MGnify metagenomic sequences, and the largest sets here
by a long way. `build-set` cuts the sources into a set under `sets/`: `fixed` is one query
set against target shards of equal size, `cross` is every query source against
every target source. Each writes a `set.tbl`, and every pipeline here takes
`--in <label>` to say which one to search.

Three benchmarks. `recall` searches every shard while sweeping nail's and
MMseqs2's prefilter sensitivity. `cloud-search` seeds once, then searches every
`(A, B)` pruning cell off those same seeds, so the pruning parameters are the
only thing moving. The grid is run twice, at `-a 5` and at `-a 0`, because a
surface at one `-a` cannot say whether a hard-pruning cell found its hits or
recovered them. `loss-decomp` asks where nail loses the hits HMMER finds,
and decomposes that across the seeding knobs: `static` against
`--mmseqs-max-seqs`, and `prog` against `--prog-n` and `--prog-f`. Every arm is
a seed list of its own and a nail that replays it; hmmer runs once outside the
sweep, because the truth set is the same for all of them and is most of the
wall clock -- 59 minutes of a 68-minute run, against 8 minutes an arm, measured
on the 2,455,939-sequence set both real labels carried before they were dropped
to a round million.

`build-set` counts each shard as it deals it and writes the counts into
`set.tbl`: counting a thousand shards afterwards is the whole deal read again,
and the count is only a metadata line. There used to be two `sizes.tbl`
formats, one per recipe, and a `build-set --in <label>` to backfill the fixed one.
All three are gone: residues are a column on the manifest like any other.

There are two table grammars, and which one a pipeline gets is settled by
whether it moves one tool's parameters in a way that changes what that tool
scores a pair. Each of the three benchmarks carries its own copy of the code
that writes and reads them, under its own `src/scores/`, holding only the half
it uses: recall has `write.rs` and `read.rs`, the two sweeps have `runs.rs`.

`recall parse scores --in <label>` reads recall into `scores.tbl`: one row per query/target
pair, one score column per tool, and a `pass` string holding one character
per run. A prefilter sweep changes which pairs a tool reports, not what it
scores them, so one column per tool is the honest shape and the cheap one --
six runs over a thousand shards is four billion rows, and a column per run is
what made an earlier grammar write 700 GB.

`cloud-search parse runs --in <label>` reads cloud-search and loss-decomp into
`runs.tbl`: one score column per run, and no pass string. A cell is the score,
or `-` where that run's seeding offered the pair and the run did not report
it, or `.` where the seeding never offered it at all. Whether a run cleared its
family's cutoff is its score against the cutoff `#= cutoffs` names, worked out
when the table is read. A seeding sweep gives every arm its own seed list, so
the offered-or-not answer differs between runs of one pair; a run's ledger row
names the seeding it replayed in a `seeds` setting, and many runs to one
seeding is the ordinary case.
`-A` and `-B` constrain the dynamic programming, so two cells can score one
pair differently, and that is what the sweep measures.
These pipelines search one shard, so the wider row costs nothing.

A run's name is stacked over its column on the dashes it is built from, so
`A10.0-B12.0-a1` takes three header lines rather than fourteen columns of
padding on every row.

Both are collected the same way: each shard on its own, written as a block in
shard order, so the file comes out the same bytes however many threads it ran.
`--threads` and `--mem` size that. `summary` asks only whether each run
cleared its family's cutoff and what the domain list says, which recall's copy
reads out of the pass string and the sweeps' out of the score columns.
`stages` reads `runs.tbl` alone, because where a pair was dropped is a
question about a run rather than about a tool.

    recall        parse scores → parse summary
    cloud-search  parse runs   → parse summary → plot
                               → parse stages
    loss-decomp   parse runs   → parse stages  → plot

`plot` draws one figure off cloud-search's summary: the `(A, B)` heatmaps of
sensitivity and wall clock, with the `full` run in the far corner. `--full-dp`
records no `-A` and no `-B`, so it has no cell; the corner is where the surface
is heading as the pruning relaxes, and it is the ceiling the pruned cells are
read against.

Each field is drawn at each `-a` and then as the difference between them, six
panels. The difference is what the second arm is run for: it says where on the
grid the sensitivity credited to `-A` and `-B` was recovered by retrying
disjoint clouds instead. A summary with one `-a` in it draws two panels.

A second figure used to scatter every cell in wall time against sensitivity.
Eighty-one points overlapping on a plane, each carrying an `-A` and a `-B` that
the position does not show, is not a figure anyone can read, so it was
deleted.

Two things sit outside the shape above, and neither asks what was found.

`cutoffs` calibrates scores. There used to be a `calibrate` that fitted what a
run would cost, over the `ladder` shape; git has it.

### The words

Naming here is deliberate and worth holding to, because the loose version of it
hides the one thing the method turns on.

- **initial reversals** -- the source sequences written backwards. These are
  not decoys. Nothing has been asked about them yet.
- **recruits** -- initial reversals that scored against some family in the
  first search. `candidate` is not the word: this repo already uses that for
  what a prefilter promotes.
- **decoys** -- recruits whose *original* has no significant match to the
  family that recruited them. Unlikely to be homologs, so they are noise, and
  noise is what a threshold is learned from.
- **rejects** -- recruits whose original *does* match. Reversed homologs
  wearing a decoy's clothes, and they are thrown out.
- **original** -- a recruit re-reversed, which is bit-identical to the source
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
Everything the calibration makes is an output of that set. What makes it a
calibration rather than a benchmark is that its product,
`data/mgy-cutoffs.tbl`, is committed and promoted by hand.

The recruits were briefly a set of their own, with a `set.tbl` and a shape.
That was wrong, and the tell was that nothing ever loaded the manifest --
`reject` globs the directory for families. `build-set` makes a set from sources
under a recipe, deterministically; what `recruit` happened to score could never
be a recipe and was never a set. Before giving anything a `set.tbl`, ask
whether `build-set` could produce it from a label.

### Rejection is a per-pair test, and that is not negotiable

**A family's null holds only sequences vetted for that family. Never pool them
across families.**

`reject` asks exactly one question about one pair: family A recruited sequence
S, so does A also match S's original? A yes makes S a reject -- the reversal of
a true member of A -- rather than a decoy for A. Family B's answer for B's own
recruit settles it for B and for nothing else.

The reason this matters more than it looks: a reversed sequence keeps a
surprisingly high similarity score against its unreversed self. That is not
only biology -- a good part of it is approximate-palindrome statistics of text,
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
cutoff would rise, and real hits would be discarded quietly -- the exact failure
rejection exists to prevent. A measurable effect from pooling is the
method reporting its own unsoundness, so treat it as a bug to find rather than
a result to keep.

The one case where B's recruits could legitimately inform A is A and B being
highly related, and that is the two families not being independent rather than
an argument for pooling.

### Pinned, for whenever cutoffs is next opened

**`--strategy fanout|union|auto`.** Two ways for `reject` to score the recruits.
`fanout` is one search per family per form, which is what the pipeline does
today. `union` is two invocations -- all of Pfam against every recruit's
reversal, then against every recruit's original -- with `learn` joining each hit
back to the family that recruited the sequence. The per-pair rule above is what
`union` must not break.

Which is faster is an open question, and **there is no usable cost model for it
yet.** nail's and mmseqs' runtimes do not scale linearly in target size, so
estimating either arm by multiplying a per-comparison rate is wrong and any
crossover derived that way is not to be trusted. The deleted `calibrate` crate
fitted the target axis with the query held fixed, which is the closest this
repository came to a cost model.

What is actually measured, and nothing more:

- `fanout` launches one process per family per direction. A per-invocation
  floor of about 0.65s was seen on a **four-family toy**, where the target was
  five sequences and contributed no measurable work. Treat it as a floor
  observed at toy scale, not a constant.
- `data/mgy-cutoffs.tbl` records what the recruits came to for MGnify: 20,795
  families, mean 5,427 recruits, max 221,046, summing to 1.13e8 (family,
  sequence) pairs. The count of *distinct* sequences those cover is not
  recorded and is the quantity `union` scales with.

`union` is clearly right at SwissProt scale, where the pool caps the union near
570k. For MGnify it is unproven either way. Do not switch MGnify over on the
strength of the SwissProt result, and do not pick between the arms from an
extrapolation -- measure both on the same input.

**`auto` needs a fitted model before it means anything**, not a rule of thumb.
Leave it unimplemented until `fanout` and `union` have been timed side by side
on real inputs.

**Measuring MGnify's union no longer waits on libsail.** 0.4.0 landed the
sequential sampler and `build-set` deals through it, so a full reversed MGnify
set is buildable rather than blocked. What it costs in practice has not been
timed.

Two things that were here are not any more, and git has both. `store import`
turned a cluster run into an ordinary run directory, reading the `.time` files
that came back with the tables, with a `rename-old-results.sh` beside it for
the older harness's filenames. `compare-scores.py` held a `scores.tbl` against
one in the older shape and said where they differed.

## benchmarks/pid

Recall as a function of the percent identity between a query and its target.
`build-set --in pid-toy|pid-real` assembles it from a profmark split: Pfam
families divided by identity, with their true targets hidden among whole TrEMBL
sequences written backwards, 100 per true pair. `truth.tbl` records which pair
is which and at what identity, and `originals.fa` holds each decoy unreversed
under the name `target.fa` gives it.

The whole benchmark is one search unit -- every query against one target file --
so the set has one row, and the truth is per pair rather than per unit. That is
why `truth` names a file instead of being a column: see the `profmark` shape.

The split itself sits at `benchmarks/pid/profmark/` rather than in the set. It
depends only on the alignments and the split parameters, both labels draw from
the same one, and the recipe points at it. Drawing it costs about a minute with
`create-profmark --onlysplit`. The toy's assembly after it took 68 s, 66 s of
that drawing 20,000 decoys from TrEMBL's 149.4 million sequences.

`pid run` searches every tool against the set, into `outputs/search/`. `pid
reject` then searches each tool against the originals of the decoys it has to
settle, into `outputs/reject/`. `parse` turns the results into the tables the
plot scripts read, and `plot` draws them. All four take `--in <label>`.

Rejection is the cutoffs rule applied per pair: a (query, decoy) hit is dropped
when the same tool, in the same mode, hits that decoy's original at E ≤ 1e-3.
A tool is searched only against the decoys one of its runs ranked at or above
that run's worst true pair. A decoy ranked below every true pair cannot move a
point of the ROC, so the curve means what it did before rejection existed.
Each tool runs at the most sensitive setting it swept. The union a tool has to
settle can be most of the decoys. On the toy with Swissprot decoys, a few of
mmseqs' true pairs scored at E ≈ 1e4, which pulled in nearly all 20,000. The
reject test compares E-values from databases of different sizes, and nothing
corrects for that yet.

`parse` has three subcommands. `recall` writes `pid.tbl`, `roc.tbl` and
`time.tbl`: one row per run per identity bin, one row per run per point of its
ROC curve, and one row per run of what it cost against what it found. `score`
correlates two named nail runs' scores, and `cells` scatters nail's cell
fraction against query by target length; both of those answer questions other
than recall against identity, and both came along with this benchmark's shape
rather than being asked of it.

`plot` names five figures -- roc, pid, time, cells, score -- and `--only`
draws a subset, one name per flag.

## benchmarks/long-seqs

How each tool's runtime scales as sequences get longer. Six paired
query/target files, where `run` searches each query against its pair and
`parse` turns that into tables; nothing draws them today. Its inputs are small and
checked in under `data/long-seqs/`, and `build-set --in long-seqs-toy` or
`--in long-seqs-real` cuts them into a set the same way every other benchmark
gets one. Its queries are sequences rather than profiles, which is what the
manifest's `query_fa` column is for, and its shape is `pairs`: one query
against one target, with the two sources separate directories so a sequence is
never on both sides.

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

# the workspace names tabl `../tabl`, which from the copy is this
ln -sfn "$(dirname "$ROOT")/tabl" "$ROOT/tmp-claude/tabl"
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
it built: `benchmarks/recall/recall run --in toy`, and so on. Every path these
crates resolve comes from a `CARGO_MANIFEST_DIR` (`util/src/tools.rs:18` for
the repo root, and pid's own for its tree), so a build in the sandbox reads the
sandbox's `data/` and
`tools/` links and writes the sandbox's `sets/`. Re-syncing an
existing copy is near instant, and the build from cold takes about fifteen
seconds.

The links are the two directories worth 4.3G between them, plus the profmark
split, which is expensive for the reasons the pid section gives. The sandbox
reads all three; nothing in it should write them, so do not run `make data` or
`make tools` from the copy. The fourth is source rather than data: `tabl` is
built from wherever that link points, so an edit to it reaches the sandbox
without a sync.

`/sets` and `/sets-toy` are excluded because a build writes them and nothing
commits them, so the copy builds its own. Without the exclude,
`--delete` removes the set a `build-set` in the sandbox just wrote, since the
real tree has nothing there to match it. long-seqs' set is checked in under
`data/long-seqs/`, which is a link, so the copy reads the real one.

`/sets` is the exclude that matters most, and dropping it is expensive rather
than wrong: a checkout that has recall's results holds about 3.5 TB there, and
rsync will copy every byte of it into the sandbox. That has happened once
already, from dropping a single exclude line, and it took the home filesystem
to 84% full. Check this list against what is on disk before running the recipe.

Edit in the real tree and re-sync, never in the sandbox: an edit in the copy is
gone at the next sync. To check an analysis against a run that finished
elsewhere, copy that run's `sets/<benchmark>/outputs/` into the sandbox and parse it
there.

## What is tracked

The downloads under `data/` and the builds under `tools/bin/` are not, and
neither is anything under `sets/` or `sets-toy/`: a set, a run and an analysis are usually
things this repository can make again.

**`sets/recall/` is the exception, and it is not backed up by being
reproducible.** It holds about 3.5 TB of results searched on a cluster -- the
hmmer set alone is 1000 shards representing some 21,926 hours of wall clock --
and nothing here can produce them again. Git ignores it like every other set,
so the only thing protecting it is that nobody deletes it. Read it,
never write it, and never point an `rsync --delete` or a `build-set --rebuild`
at that label. Its `set.tbl` was reconstructed from the
shard sizes the old harness left, and says so in its `#= imported` line: the
deal predates `set.tbl`, so the set cannot be rebuilt from its own manifest
either.

Two exceptions are committed on purpose:
`data/long-seqs/`, because nothing fetches it and ignoring it would empty a
fresh clone, and `data/mgy-cutoffs.tbl`, because learning it is a calibration
run rather than a download.

Plotting is Python and matplotlib. Three benchmarks have a `plot` command and
a `scripts/` beside their source -- cloud-search, loss-decomp and pid. recall,
cutoffs and long-seqs stop at `parse`.

Two of pid's scripts, `plot_params.py` and `plot_threads.py`, are wired to no
figure: nothing writes what they read any more.
