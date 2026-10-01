#!/usr/bin/env python3
"""Good seeds by depth in the prefilter list, family by family.

    plot_depth.py hits.tbl lists.tbl --out figures/

Two figures, one panel per target corpus each, columns square on a log rank
axis. In both, a column is the distribution across families of a fraction:
good seeds, kept over the cutoff and found by hmmer, over prefilter pairs.
Cumulative takes the family's first R pairs; binned takes only its pairs at
the ranks the column covers. A column has one cell per value the fraction can
take, up to a hundred, and a cell is coloured by how many families sit in it,
on a log scale bent so ten thousand is seven tenths of the way to all of Pfam.
An empty cell is white.
"""

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use("Agg")

import matplotlib.pyplot as plt

TOL_BLUE = "#0077BB"

# no family is white; one to all of Pfam on a log scale
BLUES = plt.get_cmap("Blues")


def table(path):
    """A `#`-headed table's rows as dicts, the header being the last `#` line
    whose first field is not a dash."""
    names, rows = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#="):
            continue
        if line.startswith("#"):
            f = line.lstrip("#").split()
            if f and not f[0].startswith("-"):
                names = f
            continue
        if line.split():
            rows.append(dict(zip(names, line.split())))
    return rows


def fractions(hits, lists, unit, run, edges, cumulative):
    """Per column edge R, each family's good-seed fraction there.

    Cumulative: good seeds among its first R pairs over R, for a family whose
    list reaches R. Binned: good seeds among its pairs from the previous edge
    to R over the pairs it has there, for a family whose list reaches past the
    previous edge. Each fraction comes with the number of pairs under it,
    which is how many values it can take."""
    length = {r["query"]: int(r["length"]) for r in lists if r["unit"] == unit and r["run"] == run}
    ranks = defaultdict(list)
    for r in hits:
        if r["unit"] == unit and r["run"] == run:
            ranks[r["query"]].append(int(r["rank"]))
    ranks = {q: np.sort(h) for q, h in ranks.items()}

    out = []
    previous = 0
    for R in edges:
        lo = 0 if cumulative else previous
        fs, ks = [], []
        for q, n in length.items():
            if n <= lo:
                continue
            pairs = min(n, R) - lo
            if cumulative and n < R:
                continue
            h = ranks.get(q)
            good = 0
            if h is not None:
                good = np.searchsorted(h, R, side="right") - np.searchsorted(h, lo, side="right")
            fs.append(good / pairs)
            ks.append(good)
        out.append((np.array(fs), np.array(ks), R - lo))
        previous = R
    return out


def draw(hits_path, lists_path, out, cumulative):
    hits = table(hits_path)
    lists = table(lists_path)
    units = list(dict.fromkeys(r["unit"] for r in lists))
    run = next(iter(dict.fromkeys(r["run"] for r in lists)))
    top = max(int(r["length"]) for r in lists if r["run"] == run)

    cells = 100
    if cumulative:
        # columns square on the log axis: edges equally spaced in log rank,
        # as many as there are 1% cells up the side
        per_decade = cells / np.log10(top)
        edges = np.unique(
            np.round(np.logspace(0, np.log10(top), int(np.log10(top) * per_decade) + 1)).astype(int)
        )
        edges = edges[edges >= 1]
    else:
        # bin widths doubling from 200, the way prog's n_take does: 1-200,
        # 201-600, 601-1,400 and so on, drawn as equal columns
        edges, width = [200], 200
        while edges[-1] < top:
            width *= 2
            edges.append(edges[-1] + width)
        edges = np.array(edges)

    fig, axes = plt.subplots(1, len(units), figsize=(8 * len(units), 6), squeeze=False,
                             layout="constrained")
    mesh = None
    if cumulative:
        # a column spans from halfway to the previous edge to halfway to the next
        xedges = np.concatenate([[edges[0]], np.sqrt(edges[:-1] * edges[1:]), [edges[-1]]])
    else:
        xedges = np.arange(len(edges) + 1)
    for ax, unit in zip(axes[0], units):
        families = len({r["query"] for r in lists if r["unit"] == unit and r["run"] == run})
        # a power scale rather than a log: a log put 50% and 100% in the
        # same shade, a linear scale put everything under 5% in one
        # log in family count, bent so that 10,000 sits at seven tenths
        # of the way to full blue rather than the 0.93 a plain log gives
        knee = min(10_000, families)
        lg = np.log10([1, knee, families])
        at = [0, 0.7, 1]
        norm = mpl.colors.FuncNorm(
            (lambda x: np.interp(np.log10(x), lg, at),
             lambda y: 10 ** np.interp(y, at, lg)),
            vmin=1, vmax=families,
        )
        columns = fractions(hits, lists, unit, run, edges, cumulative)
        if not cumulative:
            # rows are good seeds per family in the bin: 0, 1, 2, 3, then
            # ranges widening by half, up to the most any family has
            most = max((ks.max() for _, ks, _ in columns if len(ks)), default=1)
            rows = [0, 1, 2, 3, 4]
            while rows[-1] <= most:
                rows.append(int(np.ceil(rows[-1] * 1.5)))
            rows = np.array(rows)
        for c, (fs, ks, width) in enumerate(columns):
            if not len(fs):
                continue
            # over W pairs the fraction is k/W for k in 0..W, so a column
            # under 100 wide gets one cell per reachable value, centred on
            # it, rather than 1% cells most of which nothing can land in
            n = min(int(width), cells)
            if cumulative:
                ybins = np.clip((np.arange(n + 2) - 0.5) / n, 0, 1)
                counts, _ = np.histogram(fs, bins=ybins)
            else:
                counts, _ = np.histogram(ks, bins=np.append(rows, rows[-1] + 1) - 0.5)
                # equal rows, however wide the count range each one holds
                ybins = np.arange(len(rows) + 1)
            # coloured here rather than by the mesh: FuncNorm paints an
            # empty cell as if it held families, where a cell with none
            # must be nothing at all
            rgba = BLUES(norm(np.where(counts > 0, counts, 1).astype(float)))
            rgba[counts == 0, 3] = 0
            mesh = ax.pcolormesh(xedges[c : c + 2], ybins, rgba[:, None, :],
                                 edgecolors="none", antialiased=False)
            mesh.set_rasterized(True)
        if cumulative:
            ax.set_xscale("log")
            ax.set_xlim(1, top)
        else:
            ax.set_xlim(0, len(edges))
            ax.set_yticks(np.arange(len(rows)) + 0.5)
            labels = []
            for i, lo in enumerate(rows):
                hi = rows[i + 1] - 1 if i + 1 < len(rows) else None
                labels.append(f"{lo:,}" if hi is None or hi == lo else f"{lo:,}-{hi:,}")
            ax.set_yticklabels(labels, fontsize=mpl.rcParams["font.size"] * 0.8)
            ax.set_xticks(np.arange(len(edges)) + 0.5)
            lows = np.concatenate([[1], edges[:-1] + 1])
            ax.set_xticklabels([f"{lo:,}-{hi:,}" for lo, hi in zip(lows, edges)],
                               rotation=45, ha="right")
        ax.set_ylim(0, 1 if cumulative else len(rows))
        ax.set_xlabel("rank")
        ax.set_ylabel("good seeds / pairs, cumulative" if cumulative else "good seeds in the bin, per family")
        ax.set_title(unit, loc="left")
    # plain counts on the bar, with all of Pfam at the top
    ticks = [t for t in (1, 10, 100, 1000, 10000) if t < families] + [families]
    scale = mpl.cm.ScalarMappable(norm=norm, cmap=BLUES)
    bar = fig.colorbar(scale, ax=list(axes[0]), label="families", pad=0.01, ticks=ticks)
    bar.ax.set_yticklabels([f"{t:,}" for t in ticks])
    bar.ax.minorticks_off()

    name = "cumulative" if cumulative else "binned"
    path = Path(out) / f"loss-decomp-depth-{name}.pdf"
    fig.savefig(path, dpi=300)
    print(f"wrote {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("hits")
    p.add_argument("lists")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    Path(a.out).mkdir(parents=True, exist_ok=True)
    draw(a.hits, a.lists, a.out, cumulative=True)
    draw(a.hits, a.lists, a.out, cumulative=False)


if __name__ == "__main__":
    main()
