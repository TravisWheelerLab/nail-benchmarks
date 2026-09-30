#!/usr/bin/env python3
"""How far down the prefilter list the hits sit, per query, and what a static
cap pays to reach them.

    plot_depth.py depth.tbl hits.tbl --out figures/

One row per target corpus, two panels.

Left: one curve per query, the fraction of that query's hits reached by rank,
at the unbounded static seeding at the shared sensitivity. The curves flatten
at different ranks, and that spread is the case for prog seeding: a static
--mmseqs-max-seqs is one vertical line through all of them, paying for every
rank left of it under a curve that has already flattened and giving up
everything right of it under a curve still rising. The 10th, 50th and 90th
percentile curves are drawn over the web so the spread reads where the web is
solid.

Right: what a static cap costs against what it reaches. Each point is a cap at
a bin edge: the seeds nail would align, which are the seeds at or under that
rank, against the hits reached as a fraction of what hmmer found. The prog arm
is one point, at the seeds it aligned and the hits it kept.
"""

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

SCALE = 1.5
mpl.rcParams.update({"font.size": mpl.rcParams["font.size"] * SCALE})

# the palette the other benchmarks use, so figures from all of them sit together
TOL_RED = "#CC3311"
TOL_TEAL = "#009988"
TOL_BLUE = "#0077BB"
TOL_ORANGE = "#EE7733"
GREY = "#BBBBBB"

# nail's own --mmseqs-max-seqs in static mode
STATIC_DEFAULT = 300


def table(path):
    """A `#`-headed table: its `#=` lines as (key, words) and its rows as
    dicts, the header being the last `#` line whose first field is not a
    dash."""
    names, meta, rows = None, [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#="):
            f = line[2:].split()
            meta.append((f[0], f[1:]))
            continue
        if line.startswith("#"):
            f = line.lstrip("#").split()
            if f and not f[0].startswith("-"):
                names = f
            continue
        if line.split():
            rows.append(dict(zip(names, line.split())))
    return meta, rows


def per_query(hits, unit, run):
    """Each query's hit ranks, sorted, for one unit and run."""
    ranks = defaultdict(list)
    for r in hits:
        if r["unit"] == unit and r["run"] == run:
            ranks[r["query"]].append(int(r["rank"]))
    return {q: np.array(sorted(v)) for q, v in ranks.items()}


def web(ax, curves, unit):
    """One step curve per query, and the percentiles over them."""
    top = max(int(v[-1]) for v in curves.values())
    segments = []
    for ranks in curves.values():
        n = len(ranks)
        # a step at each hit, held flat to the next, and out to the edge
        x = np.concatenate([[1], np.repeat(ranks, 2), [top]])
        y = np.concatenate([[0, 0], np.repeat(np.arange(1, n + 1) / n, 2)[:-1], [1]])
        segments.append(np.column_stack([x, y]))
    ax.add_collection(
        LineCollection(segments, colors=TOL_BLUE, linewidths=0.4, alpha=0.04)
    )

    grid = np.logspace(0, np.log10(top), 300)
    reached = np.array(
        [np.searchsorted(v, grid, side="right") / len(v) for v in curves.values()]
    )
    # a low percentile of the fraction reached is a query that reaches
    # its hits late, so the 10th sits to the right of the median
    for q, style, label in [
        (10, ":", "slowest tenth of queries"),
        (50, "-", "median query"),
        (90, "--", "fastest tenth of queries"),
    ]:
        ax.plot(grid, np.percentile(reached, q, axis=0), style, color="black",
                linewidth=2, label=label)

    ax.axvline(STATIC_DEFAULT, color=TOL_RED, linewidth=1.5)
    ax.text(STATIC_DEFAULT * 1.15, 0.5, f"static default\n{STATIC_DEFAULT}",
            color=TOL_RED, fontsize=mpl.rcParams["font.size"] * 0.75)

    ax.set_xscale("log")
    ax.set_xlim(1, top)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("rank in the query's prefilter list")
    ax.set_ylabel("fraction of the query's hits reached")
    ax.set_title(f"{unit}: {len(curves):,} queries with a hit", loc="left")
    ax.legend(loc="lower right", frameon=False)


def cost(ax, depth, truth, unit, ceiling, arms):
    """Seeds aligned against sensitivity: the static caps as a curve, every
    other arm as a point."""
    caps = [r for r in depth if r["unit"] == unit and r["run"] == ceiling and r["lo"] != "-"]
    seeds = np.cumsum([int(r["seeds"]) for r in caps])
    hits = np.cumsum([int(r["hits"]) for r in caps]) / truth
    ax.plot(seeds, hits, "-o", color=TOL_RED, markersize=6, linewidth=2,
            label=f"static cap at {ceiling}")
    # a cap's rank beside its point, skipped where the points crowd
    last_x = 0
    for r, x, y in zip(caps, seeds, hits):
        if x < last_x * 1.04:
            continue
        last_x = x
        ax.annotate(r["hi"], (x, y), textcoords="offset points", xytext=(6, -12),
                    fontsize=mpl.rcParams["font.size"] * 0.6, color=TOL_RED)

    for run, color, marker in arms:
        rows = [r for r in depth if r["unit"] == unit and r["run"] == run and r["lo"] != "-"]
        if not rows:
            continue
        x = sum(int(r["seeds"]) for r in rows)
        y = sum(int(r["hits"]) for r in rows) / truth
        ax.scatter([x], [y], s=110, color=color, marker=marker, zorder=3, label=run)

    ax.set_xscale("log")
    ax.set_xlabel("seeds nail aligns")
    ax.set_ylabel("hits reached / hmmer's hits")
    ax.set_title(f"{unit}: what a cap pays", loc="left")
    ax.legend(loc="lower right", frameon=False)


def draw(depth_path, hits_path, out):
    meta, depth = table(depth_path)
    _, hits = table(hits_path)

    truth = {words[0]: int(words[1]) for key, words in meta if key == "truth"}
    units = list(dict.fromkeys(r["unit"] for r in depth))
    runs = list(dict.fromkeys(r["run"] for r in depth))

    # the ceiling is the unbounded static arm at the shared sensitivity; the
    # rest are drawn as points against it
    static = [r for r in runs if r.startswith("s")]
    ceiling = static[0]
    others = [(r, TOL_ORANGE, "s") for r in static[1:]]
    others += [(r, TOL_TEAL, "D") for r in runs if r.startswith("prog")]

    fig, axes = plt.subplots(len(units), 2, figsize=(16, 6.5 * len(units)), squeeze=False)
    for (left, right), unit in zip(axes, units):
        web(left, per_query(hits, unit, ceiling), unit)
        cost(right, depth, truth[unit], unit, ceiling, others)

    fig.tight_layout()
    path = Path(out) / "loss-decomp-depth.pdf"
    fig.savefig(path)
    print(f"wrote {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("depth")
    p.add_argument("hits")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    Path(a.out).mkdir(parents=True, exist_ok=True)
    draw(a.depth, a.hits, a.out)


if __name__ == "__main__":
    main()
