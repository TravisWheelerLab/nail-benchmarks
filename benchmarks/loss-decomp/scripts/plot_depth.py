#!/usr/bin/env python3
"""How good the seeds are at each depth of the prefilter list, family by
family.

    plot_depth.py hits.tbl lists.tbl --out figures/

One panel per target corpus. The x axis is depth in a query's prefilter list,
in bins that double from the first; the y axis is a family's yield in that
bin: of its prefilter pairs at those ranks, the fraction that became hits, kept
over the cutoff and found by hmmer. Every family whose list reaches the bin is
one sample, so the column over a bin is the distribution of yield across
families, drawn as a density on a log axis with the median and the 10th and
90th percentiles over it.

The yield means the same thing for a family with 3 hits and one with 3,000,
which is what lets the families be pooled. At rank 1 it is 0 or 1 per family
and says nothing; by rank 200 it is a band, and how that band falls with depth,
and how wide it is at each depth, is what a stopping rule is up against: a
static cap stops every family at one depth, prog stops each where its own
yield gives out.

Families with a yield of exactly zero in a bin sit off a log axis, so their
share is written under each column instead.
"""

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

SCALE = 1.5
mpl.rcParams.update({"font.size": mpl.rcParams["font.size"] * SCALE})

# the palette the other benchmarks use, so figures from all of them sit together
TOL_RED = "#CC3311"
TOL_BLUE = "#0077BB"

# nail's own --mmseqs-max-seqs in static mode
STATIC_DEFAULT = 300

# yields below this are drawn at the axis floor; the strip under it holds
# the share of families at exactly zero
FLOOR = 1e-5
STRIP = FLOOR / 5


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


def edges(width, top):
    """The bins: (lo, hi) doubling from `width` until `top` is covered."""
    out, lo, hi = [], 1, width
    while lo <= top:
        out.append((lo, hi))
        lo, hi = hi + 1, hi * 2
    return out


def yields(hits, lists, unit, run, bins):
    """Per bin, each family's yield there: hits over prefilter pairs, for
    every family whose list reaches the bin."""
    length = {r["query"]: int(r["length"]) for r in lists if r["unit"] == unit and r["run"] == run}
    hit = defaultdict(lambda: np.zeros(len(bins), dtype=int))
    for r in hits:
        if r["unit"] != unit or r["run"] != run:
            continue
        rank = int(r["rank"])
        for b, (lo, hi) in enumerate(bins):
            if lo <= rank <= hi:
                hit[r["query"]][b] += 1
                break

    out = []
    for b, (lo, hi) in enumerate(bins):
        ys = []
        for q, n in length.items():
            if n < lo:
                continue
            pairs = min(n, hi) - lo + 1
            ys.append(hit[q][b] / pairs)
        out.append(np.array(ys))
    return out


def panel(ax, per_bin, bins, unit):
    """One column of density per bin, the percentiles over it, and the share
    of families at zero under it."""
    ybins = np.logspace(np.log10(FLOOR), 0, 26)
    grid = np.zeros((len(ybins) - 1, len(bins)))
    for b, ys in enumerate(per_bin):
        nonzero = ys[ys > 0]
        if len(nonzero):
            grid[:, b], _ = np.histogram(np.clip(nonzero, FLOOR, 1), bins=ybins)

    x = np.arange(len(bins) + 1)
    mesh = ax.pcolormesh(x, ybins, np.where(grid > 0, grid, np.nan),
                         cmap="Blues", norm=LogNorm(vmin=1, vmax=max(grid.max(), 2)))

    # percentiles over the families with a hit in the bin, which is what
    # the density above draws; the rest are the share written under it
    centers = x[:-1] + 0.5
    for q, style, label in [
        (90, "--", "90th percentile family with a hit"),
        (50, "-", "median family with a hit"),
        (10, ":", "10th percentile family with a hit"),
    ]:
        v = np.array([
            np.percentile(ys[ys > 0], q) if np.any(ys > 0) else np.nan for ys in per_bin
        ])
        ax.plot(centers, v, style, color="black", linewidth=2, label=label)

    # the share of families with no hit at all in the bin, in the strip
    # under the axis floor
    ax.axhline(FLOOR, color="black", linewidth=0.8)
    for b, ys in enumerate(per_bin):
        if len(ys) == 0:
            continue
        zero = np.mean(ys == 0)
        ax.text(centers[b], STRIP * 1.15, f"{zero:.0%}\nat 0\nof {len(ys):,}",
                ha="center", va="bottom", fontsize=mpl.rcParams["font.size"] * 0.6)

    # where nail's static default would stop every family
    for b, (lo, hi) in enumerate(bins):
        if lo <= STATIC_DEFAULT <= hi:
            at = b + (np.log2(STATIC_DEFAULT / lo + 1) if b else STATIC_DEFAULT / hi)
            ax.axvline(at, color=TOL_RED, linewidth=1.5)
            ax.text(at + 0.1, 0.5, f"static default {STATIC_DEFAULT}", color=TOL_RED,
                    fontsize=mpl.rcParams["font.size"] * 0.7)
            break

    ax.set_yscale("log")
    ax.set_ylim(STRIP, 1.05)
    ax.set_xlim(0, len(bins))
    ax.set_xticks(centers)
    ax.set_xticklabels([f"{lo:,}\n{hi:,}" for lo, hi in bins], fontsize=mpl.rcParams["font.size"] * 0.6)
    ax.set_xlabel("depth in the query's prefilter list")
    ax.set_ylabel("yield: hits / prefilter pairs at that depth, per family")
    ax.set_title(unit, loc="left")
    ax.legend(loc="upper right", frameon=False, fontsize=mpl.rcParams["font.size"] * 0.8)
    return mesh


def draw(hits_path, lists_path, out):
    meta, hits = table(hits_path)
    _, lists = table(lists_path)

    units = list(dict.fromkeys(r["unit"] for r in lists))
    runs = list(dict.fromkeys(r["run"] for r in lists))
    # the unbounded static arm at the shared sensitivity, the first arm
    ceiling = runs[0]
    width = 200
    top = max(int(r["length"]) for r in lists if r["run"] == ceiling)
    bins = edges(width, top)

    fig, axes = plt.subplots(1, len(units), figsize=(9.5 * len(units), 8), squeeze=False,
                             layout="constrained")
    mesh = None
    for ax, unit in zip(axes[0], units):
        mesh = panel(ax, yields(hits, lists, unit, ceiling, bins), bins, unit)
    fig.colorbar(mesh, ax=list(axes[0]), label="families", pad=0.01, shrink=0.8)

    fig.suptitle(f"{ceiling}, static, unbounded: how good the seeds are by depth", x=0.01, ha="left")
    path = Path(out) / "loss-decomp-depth.pdf"
    fig.savefig(path)
    print(f"wrote {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("hits")
    p.add_argument("lists")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    Path(a.out).mkdir(parents=True, exist_ok=True)
    draw(a.hits, a.lists, a.out)


if __name__ == "__main__":
    main()
