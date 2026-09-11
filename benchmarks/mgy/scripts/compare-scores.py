#!/usr/bin/env python3
"""Hold a scores.tbl written in the new shape against one written in the old.

The two files say the same things differently: the old one carries a score and
a pair of cutoff columns per run, the new one a score per tool and a pass
character per run. Both are read down to the same three sets --

    (query, target, run, passed)
    (query, target, tool, score)
    (query, target, sorted domain scores)

-- which have to be identical. The new file is also checked for the order it
promises: sorted inside a block, one block per shard, shards as the `#= target`
lines list them.

    compare-scores.py old.tbl new.tbl [-c 2]

`-c` is the cutoff column the old file was written with, which it does not
record; the new one does, and a mismatch is reported.
"""

import argparse
import sys
from collections import defaultdict

TOOLS = ("nail", "mmseqs", "hmmer")

# the old table holds hmmer to nail's cutoff, as the calibration does
CUTOFF = {"nail": "cut_nail", "hmmer": "cut_nail", "mmseqs": "cut_mmseqs"}


def runs_of(lines):
    """The `#= run` lines: a name and a tool each, in order."""
    out = []
    for line in lines:
        if line.startswith("#= run "):
            fields = line.split()
            out.append((fields[2], fields[3]))
    return out


def read_old(path):
    """The old shape: one score column per run, two cutoff columns per row."""
    text = open(path).read().splitlines()
    runs = runs_of(text)

    header = next(
        line[1:].split()
        for line in text
        if line.startswith("#") and not line.startswith("#=") and " query " in f" {line[1:]} "
    )

    passed, scores, doms = {}, {}, {}
    disagree = []

    for line in text:
        if line.startswith("#") or not line.strip():
            continue

        cells = dict(zip(header, line.split()))
        pair = (cells["query"], cells["target"])

        for name, tool in runs:
            score = cells[name]
            cut = cells[CUTOFF[tool]]

            passed[(pair, name)] = (
                score != "-" and cut != "-" and float(score) >= float(cut)
            )

            if score == "-":
                continue

            # every run of one tool should have given the pair the same score,
            # which is the claim the new table's single column rests on
            seen = scores.setdefault((pair, tool), score)
            if seen != score:
                disagree.append((pair, tool, seen, score))

        # one `<run>_dom` column per hmmer run, comma-joined and best first
        tail = [cells[f"{name}_dom"] for name, tool in runs if tool == "hmmer"]
        doms[pair] = sorted(
            float(x) for cell in tail if cell != "-" for x in cell.split(",")
        )

    return runs, passed, scores, doms, disagree


def read_new(path):
    """The new shape: a pass character per run, a score column per tool."""
    text = open(path).read().splitlines()

    if not text or text[0] != "#= format scores 2":
        sys.exit(f"{path} does not open `#= format scores 2`")

    runs = runs_of(text)
    legend = next((line.split()[2:] for line in text if line.startswith("#= pass ")), [])
    shards = [line.split()[2] for line in text if line.startswith("#= target ")]

    header = next(
        line[1:].split()
        for line in text
        if line.startswith("# query")
    )
    tools = [h for h in header if h in TOOLS]

    passed, scores, doms = {}, {}, {}
    order, seen, block, rows = [], set(), None, 0

    for line in text:
        if line.startswith("#= shard "):
            block = line.split()[2]
            if block in seen:
                sys.exit(f"{path} opens shard {block} twice")
            seen.add(block)
            order.append(block)
            last = None
            continue

        if line.startswith("#= end "):
            if int(line.split()[2]) != rows:
                sys.exit(f"{path} says {line.split()[2]} rows and holds {rows}")
            continue

        if line.startswith("#") or not line.strip():
            continue

        rows += 1
        cells = line.split()
        pair = (cells[0], cells[1])

        if last is not None and pair <= last:
            sys.exit(f"{path} is out of order at {pair} in shard {block}")
        last = pair

        for at, (name, _) in enumerate(runs):
            passed[(pair, name)] = cells[2][at].isupper()

        for at, tool in enumerate(tools):
            if cells[3 + at] != "-":
                scores[(pair, tool)] = cells[3 + at]

        tail = cells[4 + len(tools)]
        doms[pair] = sorted(float(x) for x in tail.split(",")) if tail != "-" else []

    if order != shards:
        sys.exit(f"{path} blocks are {order}, the `#= target` lines are {shards}")

    if legend != [name for name, _ in runs]:
        sys.exit(f"{path} legend is {legend}, the runs are {[n for n, _ in runs]}")

    return runs, passed, scores, doms


def differences(what, old, new, most=5):
    keys = set(old) | set(new)
    bad = [k for k in keys if old.get(k) != new.get(k)]

    if not bad:
        print(f"{what}: {len(keys)} identical")
        return 0

    print(f"{what}: {len(bad)} of {len(keys)} differ")
    for key in sorted(bad, key=str)[:most]:
        print(f"  {key}: old {old.get(key)!r}, new {new.get(key)!r}")

    return len(bad)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("old")
    parser.add_argument("new")
    parser.add_argument("-c", type=int, default=2)
    args = parser.parse_args()

    old_runs, old_passed, old_scores, old_doms, disagree = read_old(args.old)
    new_runs, new_passed, new_scores, new_doms = read_new(args.new)

    if old_runs != new_runs:
        sys.exit(f"the files declare different runs:\n  {old_runs}\n  {new_runs}")

    if disagree:
        print(f"note: {len(disagree)} pair(s) scored differently by two runs of one tool")
        for pair, tool, a, b in disagree[:5]:
            print(f"  {pair} {tool}: {a} and {b}")

    pairs = differences(
        "pairs",
        {p: True for p, _ in old_passed},
        {p: True for p, _ in new_passed},
    )

    bad = pairs
    bad += differences("pass flags", old_passed, new_passed)
    bad += differences("tool scores", old_scores, new_scores)
    bad += differences("domain scores", old_doms, new_doms)

    print("identical" if bad == 0 else f"{bad} differences")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
