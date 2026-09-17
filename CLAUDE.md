# nail-benchmarks

Benchmarks for [nail](https://github.com/travisWheelerLab/nail), a profile HMM
search tool, against the tools it is compared to: HMMER, MMseqs2, BLAST, LAST
and DIAMOND.

A Cargo workspace at the repo root. Binaries are the interface. There is one
per benchmark, and what each of them builds, runs and works out goes in the
store.

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
store/                    what a build, a run or an analysis produced
benchmarks/shim           the build-and-run shim every binary links to

benchmarks/util/          set, store, ledger, tbl, tools, split, cut, ...
benchmarks/search/        the tool command builders every pipeline composes
benchmarks/scores/        the pair tables, the analyses, and parse

benchmarks/build-set/     cuts sources into a set: fixed | ladder | pairs
benchmarks/store/         what can be done to the store: import | clean

benchmarks/recall/        recall against prefilter sensitivity   [fixed]
benchmarks/cloud-search/  the (A, B) pruning surface             [fixed]
benchmarks/loss-decomp/   where nail loses hmmer's hits, by stage [fixed]
benchmarks/calibrate/     what a run will cost  [ladder]  ON ICE, see below
benchmarks/cutoffs/       per-family score cutoffs from decoys    [fixed]
benchmarks/long-seqs/     how the tools scale with length         [pairs]
benchmarks/pid/           recall against percent identity          [profmark]
```

One binary per benchmark, and the shape in brackets is the set it reads. Three
libraries sit under them, and two binaries that belong to no benchmark:
`build-set` makes the sets, `store` handles what a run left behind.

`benchmarks/mgy/` is empty. It held the `inputs/` and `outputs/` trees the old
`mgy` crate wrote before the store existed, and those moved into the store on
16 September 2026 -- see the warning under **What is tracked**, which is the
one thing in this file to read before running anything near `store/`.

### How much of the libraries is actually shared

Counted per public item, by how many binaries reference it. A part with one
consumer belongs in that consumer rather than in a library, and the count is
the test rather than how well the library reads from the inside.

About a sixth of the two libraries has a single consumer, and all of it is
recall's:

- `search`, 428 lines. Roughly 87 are recall alone -- `NAIL_S`, `MMSEQS_S`,
  `MMSEQS_MAX_SEQS`, `createdb` and the whole `Mmseqs` builder, since recall is
  the only benchmark that sweeps mmseqs. Another 16 have no consumer at all:
  `cat`, `MMSEQS_K`, `Dirs.results`. The remaining ~325 -- `Dirs`, `Split`,
  `Hmmer`, `Bins`, `jobs`, `HMMER_CPU` -- are used by recall, cloud-search and
  loss-decomp alike.
- `scores`, 4,200 lines. Roughly 680 are recall alone: `write.rs` and
  `read.rs`, which are the `scores.tbl` grammar. Nothing at all is exclusive to
  cloud-search or to loss-decomp.

The seam that shows up is the one the two table grammars already draw, rather
than a split between reading and analysing: `scores.tbl` is recall's and
nobody else's, `runs.tbl` is cloud-search's and loss-decomp's. loss-decomp's use of
the crate is a strict subset of cloud-search's -- both write `runs.tbl` and
read `stages`, and cloud-search also reads `summary`.

All three reach `scores` only through `scores::parse`. No binary names
`analyze`, `write`, `runs`, `frame`, `read`, `shard` or `collect`, so what sits
behind that one door can move without touching a caller.

`cutoffs` is the benchmark that uses neither library. It runs the same three
tools and builds every command by hand, at `cutoffs/src/main.rs:471` and after.

No benchmark looks on `PATH`. `util::tools` holds the path to every binary and
every download, and a benchmark reads it rather than guessing, so a run uses
whichever `hmmsearch` was built for this repo.

`pid` reads a set like everything else now. What kept it out was a
record-level truth table, and the answer was for `set.tbl` to name the file
rather than carry it -- see the `profmark` shape. It still keeps its own copy of
the command builders in `src/search.rs`, and its `profmark/` split at the crate
root.

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
  tables nail, HMMER, MMseqs2 and BLAST produce.
- `tabl` writes the padded, `#`-headed tables. It is not published: the
  workspace takes it as a path dependency on a sibling checkout, `../tabl`.
- `feisty` sits in `[workspace.dependencies]` and no member uses it.

## The shape a benchmark has

Three stages, joined by data rather than by code, and one tree holding what
each of them produced. Where an artifact lives follows from what it is, not
from which crate made it, so a set built by one benchmark is readable by
another without either naming the other's directory.

A set and everything derived from it share a directory, so what a run was
searched against is the directory it sits in:

```
store/<set>/
├── inputs/
│   ├── set.tbl                  one row per search unit, and what it is made of
│   └── ...                      the queries and targets it names
├── outputs/<run>/
│   ├── manifest.tbl             every command, its wall clock, its exit code
│   ├── ledger.tbl               one row per run per shard, and what it cost
│   └── results/<run>.<shard>.tbl   one hit table per run, per target shard
├── analysis/<run>/              what parse worked out, and the figures
└── tmp/<run>/                   scratch, and nothing worth keeping
```

Nothing in the code knows that layout. It is what the `paths.toml` files happen
to say, and moving a set to a scratch disk is editing a line rather than
changing anything that compiles.

The scratch sits beside the record rather than inside it, so
`outputs/<run>/` holds the record and only the record, and the whole of `tmp/`
can go at any time without touching it. The search pipelines still
take `--tmp` to put their own scratch somewhere else, a scratch disk being the
usual reason; `build` and `cutoffs` have no such flag and never had one.

`store clean --paths <crate>/paths.toml --in <label>` removes what that label
names as produced -- the run, the analysis and the scratch -- after printing how
many files and how many bytes are about to go and waiting for a `y`. The set
goes only with `--all`, since a build is expensive. Never point it at recall's
`real` label: see **What is tracked**.

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
read a set's shape the same way. `fixed` and `ladder` produce the same table
with different attribute columns, and a pipeline reading it does not learn which
one ran.

Both are written by the producer in the pass that produces the artifact, and
nothing else ever writes them. That discipline is the whole defence against a
manifest that says one thing while the directory says another.

### The shapes

An open manifest lets one format describe a deal of shards and a nest of rungs.
The cost is that a set and the benchmark reading it can disagree with nothing
saying so, and a `ladder` searched as if it were `fixed` is a sweep over the
product of both axes reported as a list of targets. A shape closes that: the
builder stamps `#= shape` and the benchmark names the one it reads, and
`Set::load_as` holds them to each other before a tool runs.

| shape | representations | attributes | read by |
|---|---|---|---|
| `fixed` | `query_hmm`, `query_sto`, `query_db`, `target` | `shard`, `seqs`, `residues`, `bytes` | recall, cloud-search, loss-decomp |
| `reversed` | the same as `fixed` | the same as `fixed` | cutoffs |
| `ladder` | the same | `query_rung`, `target_rung`, `query_residues`, `target_residues` | calibrate |
| `pairs` | `query_fa`, `target` | `pair`, `query_residues`, `residues` | long-seqs |
| `profmark` | `query_hmm`, `query_sto`, `query_fa`, `target` | `truth` | pid |

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
columns name theirs. The shape check sees that the file is named, not what is
in it.

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

## paths.toml

No crate resolves a location. A tool is told where its inputs are and where its
outputs go, and the telling is a `paths.toml` in the tool's own crate directory,
one table per label:

```toml
[toy]
set      = "../../store/toy/inputs"
run      = "../../store/toy/outputs/recall"
analysis = "../../store/toy/analysis/recall"
tmp      = "../../store/toy/tmp/recall"
```

A label is a whole set of paths under one name, so a toy run and a real run
differ by a word: `recall run --in toy`. Running a tool without `--in` prints
the labels its file holds.

Every benchmark has exactly two, `toy` and `real`, over two sets of its own:
`build-set` names them `<benchmark>-toy` and `<benchmark>-real`. A benchmark
sharing a set with another needs a reason, because a shared set hides what a
benchmark reads. cloud-search and loss-decomp search a single unit, and they did
that by opening shard 1 of recall's thousand, so their real sets are one shard
the size of one of recall's rather than a thousand of them.

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
since nothing here ran those commands and an exit code or an argv would have to
be invented. They arrive with `.time` files beside the tables, and
`store import` reads those. It is the only command that writes a ledger by
hand.

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
- `nail` reads the one column of nail's table that `libsail`'s layout does not
  carry.
- `time` reads the `.time` file of a run made outside the harness, in any of
  the six formats a `time` command might have written.

## The Pfam-against-MGnify benchmarks

Pfam profiles against MGnify metagenomic sequences, and the largest sets here
by a long way. `build-set` cuts the sources into a set under `store/`:
`fixed` is one query set against target shards of equal size, `ladder` is
nested rungs on both axes, each a prefix of the one above. Each writes a
`set.tbl`, and every pipeline here takes `--in <label>` to say which one to
search.

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
scores a pair.

`recall parse scores --in <label>` reads recall into `scores.tbl`: one row per query/target
pair, one score column per tool, and a `pass` string holding one character
per run. A prefilter sweep changes which pairs a tool reports, not what it
scores them, so one column per tool is the honest shape and the cheap one --
six runs over a thousand shards is four billion rows, and a column per run is
what made an earlier grammar write 700 GB.

`cloud-search parse runs --in <label>` reads cloud-search and loss-decomp into `runs.tbl`: one score
column per run, plus a `seeded` column where the pipeline kept a seed list.
`-A` and `-B` constrain the dynamic programming, so two cells can score one
pair differently, and that is what the sweep measures.
These pipelines search one shard, so the wider row costs nothing.

Both are collected the same way: each shard on its own, written as a block in
shard order, so the file comes out the same bytes however many threads it ran.
`--threads` and `--mem` size that. `summary` is one streaming pass and reads
either grammar -- it reads only the `pass` string and the domain list, and
those sit in the same place in both. `stages` reads `runs.tbl` alone, because where a
pair was dropped is a question about a run rather than about a tool.

    recall        parse scores → parse summary
    cloud-search  parse runs   → parse summary → plot
    loss-decomp   parse runs   → parse stages

`plot` draws one figure off cloud-search's summary: the `(A, B)` heatmaps of
sensitivity and wall clock, with the `full` run in the far corner. `--full-dp`
records no `-A` and no `-B`, so it has no cell; the corner is where the surface
is heading as the pruning relaxes, and it is the ceiling the pruned cells are
read against.

A second figure used to scatter every cell in wall time against sensitivity.
Eighty-one points overlapping on a plane, each carrying an `-A` and a `-B` that
the position does not show, is not a figure anyone can read, so it was
deleted.

Three things sit outside the shape above, and none of them asks what was
found.

**`calibrate` is on ice: leave it alone until Jack says otherwise.** It may
well be deleted. Do not fix, tidy, extend or test it, and do not count it when
working out how many benchmarks use a shared crate -- a symbol only calibrate
and one other crate reference has one consumer, not two. The rest of this
section describes it as it stands.

`calibrate` asks what a run will cost. `calibrate run` times the searches
the benchmarks are built out of -- nail's seeding, the alignment off those
seeds, `mmseqs search` and `hmmsearch` -- over the target rungs of the ladder,
building every command through the same `search` helpers the benchmarks use,
so what is timed is what will run. `--parts` narrows it to a subset of those
four. `calibrate fit` turns the timings into `cost.tbl`; `calibrate predict`
composes a pipeline out of them, folding the way the ledger does so the total
is wall clock rather than core-seconds.

Only the searches are timed. Cutting the query up, building mmseqs' database
and reformatting what it found are real wall clock, and none of them is what
the benchmarks compare.

The query is every Pfam family at every rung, because that is what the
benchmarks search with. So the only thing that moves is the target, and a
search costs

```text
intercept + slope * target_residues
```

An earlier version swept a grid on both axes and fitted a query term, a target
term and their product. Holding the query fixed collapses three of those into
the two here, and the terms it drops were the ones carrying the error: the
product coefficient did most of the work at the sizes being predicted and was
the least determined thing in the model, so dropping a single rung from the
grid moved nail's predicted cost by a factor of two.

What is left is an intercept that is mostly the query and a slope that is
entirely the target. The intercept is large, since nail builds an mmseqs
profile database out of 20,795 HMMs on every invocation, and at the bottom of
the ladder it is nearly the whole cost: over the first three rungs a fourfold
increase in target size did not move the total past the run-to-run noise. The
slope becomes measurable only once the target term clears that noise, and that
is why the rungs double all the way to 128,000 sequences. The rung of a single
sequence at the bottom measures the intercept on its own.

A whole nail search is the seeding plus the alignment off those seeds, so the
end-to-end run is not timed. Over a ladder spanning 131x the two halves came
to within 1.1% of it at every sensitivity, so the third timing was dropped.

Every `mmseqs search` here passes `-k 6`, which is what nail passes the
prefilter it seeds with. mmseqs' own default is 0, meaning it picks a k-mer
length from the size of the database, so a run left on the default would
search with a different k at every rung and a different one again from the k
inside nail.

`fit` holds the top rung back, fits on the rest, and scores its own prediction
of the rung it did not see, so a model that extrapolates badly reports it in a
`holdout` column instead of being believed. That column is the one to read
first. On the eight-rung ladder every part over-predicted, from 2.4% for the
alignments to 25% for mmseqs at its lower sensitivity, so a prediction off it
reads as an upper bound. The alignments are the rows to trust: their cost is
almost all target work, and they are the only ones the ladder pinned to better
than 3%.

`cutoffs` is the other calibration, and the words do not mean the same thing:
cutoffs calibrates scores, calibrate calibrates cost.

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
crossover derived that way is not to be trusted. `calibrate` is the machinery
for fitting this properly, and even it only fits the target axis with the query
held fixed.

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

**Measuring MGnify's union need not wait.** Recruiting a single reversed shard
and counting distinct recruits gives an estimate without building a full
reversed set, which is otherwise blocked on libsail 0.4.0's sequential sampler.
That sampler is the streaming `Selector` sketched in
`tmp-claude/libsail-sampling-proposal.md` rather than the `Indexable::sample`
that proposal leads with; combined with buffered per-shard output it should
build a real MGnify set orders of magnitude faster than permuting 2.4 billion
indices.

`store import` writes a `ledger.tbl` for result tables produced elsewhere, out of
the `.time` files that came back with them, which turns a search run on a
cluster into an ordinary run directory in the store. It takes `--set` so it can
check the shards it was handed against the ones the build made.
`benchmarks/store/scripts/rename-old-results.sh` renames the older harness's
files into the names it expects.

One dev tool sits beside all of that and belongs to no pipeline.
`benchmarks/scores/scripts/compare-scores.py` reduces a `scores.tbl` and one in
the older shape to the same sets of pairs, tool scores and pass flags, and says
where they differ.

## benchmarks/pid

Recall as a function of the percent identity between a query and its target.
`build-set --in pid-toy|pid-real` assembles it from a profmark split: Pfam
families divided by identity, their true targets hidden in a Swissprot decoy
background, and `truth.tbl` recording which pair is which and at what identity.

The whole benchmark is one search unit -- every query against one target file --
so the set has one row, and the truth is per pair rather than per unit. That is
why `truth` names a file instead of being a column: see the `profmark` shape.

The split itself sits at `benchmarks/pid/profmark/` rather than in the set. It
depends only on the alignments and the split parameters, both labels draw from
the same one, and the recipe points at it. Drawing it costs about a minute with
`create-profmark --onlysplit`; the assembly after it is seconds.

`pid run` searches every tool against the set, `parse` turns the results into
the tables the plot scripts read, and `plot` draws them. All three take
`--in <label>`.

## benchmarks/long-seqs

How each tool's runtime scales as sequences get longer. Six paired
query/target files, where `run` searches each query against its pair and
`parse` turns that into the plot scripts' tables. Its inputs are small and
checked in under `data/long-seqs/`, and `build-set --in long-seqs-toy` or
`--in long-seqs-real` cuts them into a set the same way every other benchmark
gets one. Its queries are sequences rather than profiles, which is what the
manifest's `query_fa` column is for, and its shape is `pairs`: one query
against one target, with the two sources separate directories so a sequence is
never on both sides.

## Testing behaviour: work in your own copy

This working tree belongs to whoever is at the keyboard. A pipeline writes into
`store/`, a `build` subcommand writes a set there, and when two people write
there at once neither can tell which results are theirs. So run nothing here.
Copy the source into `tmp-claude/sandbox/`, link the expensive directories, and
work in the copy.

```bash
ROOT=$(git rev-parse --show-toplevel)
SB=$ROOT/tmp-claude/sandbox

rsync -a --delete \
  --exclude '.git/' --exclude 'target/' --exclude 'tmp-claude/' \
  --exclude '/data' --exclude '/tools' --exclude '/benchmarks/pid/profmark' \
  --exclude '/store' \
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
touching the sandbox's own `store/`, `target/`, pid `inputs/` or links, since
rsync leaves excluded paths on the receiving side alone. The link paths and
`/store` are excluded without a trailing slash on purpose: a pattern ending in
`/` matches only a directory, and on the receiving side three of these are
symlinks, so `--delete` removes them. It copies uncommitted edits, which is the
point: what wants testing is usually not committed yet.

The sandbox gets its own `store/`, which is what makes a build or a run in the
copy safe: `util::tools::repo()` resolves from `util`'s own
`CARGO_MANIFEST_DIR`, so a build inside the copy resolves the copy's root and
writes the copy's store.

Then run inside `$SB`, through its own shims. Every binary links to the same
`benchmarks/shim`, which builds the crate the link is named after and runs what
it built: `benchmarks/recall/recall run --in toy`, and so on. Every path these
crates resolve comes from a `CARGO_MANIFEST_DIR` (`util/src/tools.rs:18` for
the repo root, and pid's own for its tree), so a build in the sandbox reads the
sandbox's `data/` and
`tools/` links and writes the sandbox's `store/`. Re-syncing an
existing copy is near instant, and the build from cold takes about fifteen
seconds.

The links are the two directories worth 4.3G between them, plus the profmark
split, which is expensive for the reasons the pid section gives. The sandbox
reads all three; nothing in it should write them, so do not run `make data` or
`make tools` from the copy. The fourth is source rather than data: `tabl` is
built from wherever that link points, so an edit to it reaches the sandbox
without a sync.

`/store` and pid's `inputs/` are excluded for the same reason: a build writes
them and nothing commits them, so the copy builds its own. Without the exclude,
`--delete` removes the set a `build-set` in the sandbox just wrote, since the
real tree has nothing there to match it. long-seqs' set is checked in under
`data/long-seqs/`, which is a link, so the copy reads the real one.

`/store` is the exclude that matters most, and dropping it is expensive rather
than wrong: a checkout that has the mgy results holds about 3.5 TB there, and
rsync will copy every byte of it into the sandbox. That has happened once
already, from dropping a single exclude line, and it took the home filesystem
to 84% full. Check this list against what is on disk before running the recipe.

`/benchmarks/mgy/inputs/` is still excluded and now costs nothing, since that
tree is empty. It is kept so that a checkout made before the move syncs the
same way.

Edit in the real tree and re-sync, never in the sandbox: an edit in the copy is
gone at the next sync. To check an analysis against a run that finished
elsewhere, copy that run's `store/runs/<run>/` into the sandbox and parse it
there.

## What is tracked

The downloads under `data/` and the builds under `tools/bin/` are not, and
neither is anything under `store/`: a set, a run and an analysis are usually
things this repository can make again.

**`store/recall-real/` is the exception, and it is not backed up by being
reproducible.** It holds about 3.5 TB of results searched on a cluster -- the
hmmer set alone is 1000 shards representing some 21,926 hours of wall clock --
and nothing here can produce them again. It is ignored by git like the rest of
the store, so the only thing protecting it is that nobody deletes it. Read it,
never write it, and never point `store clean`, an `rsync --delete` or a
`build-set --rebuild` at that label. Its `set.tbl` was reconstructed from the
shard sizes the old harness left, and says so in its `#= imported` line: the
deal predates `set.tbl`, so the set cannot be rebuilt from its own manifest
either.

Two exceptions are committed on purpose:
`data/long-seqs/`, because nothing fetches it and ignoring it would empty a
fresh clone, and `data/mgy-cutoffs.tbl`, because learning it is a calibration
run rather than a download.

Plotting is Python and matplotlib, under each benchmark's `scripts/`.
