#!/usr/bin/env python3
"""The cloud-search figures, out of the summary.tbl `mgy parse summary` writes.

    plot.py summary.tbl --out figures/

One row per run, and the settings a run was swept over are columns. These want
-A and -B among them, so they draw for cloud-search and skip themselves for the
pipelines that swept something else.

One figure. `heatmaps` is the surface -- sensitivity and wall time over the
(A, B) grid -- with the unpruned `--full-dp` run in the far corner. It records
no -A and no -B, so it has no cell of its own; the corner is where the surface
is heading as the pruning relaxes, and it is the ceiling every pruned cell is
read against.

A table holding two -a draws each field at each of them and then the difference,
six panels in two rows. The difference panel is the one to read: it says where
on the grid the sensitivity credited to -A and -B was recovered by retrying
disjoint clouds instead. A table holding one -a draws the two panels alone, so
a summary written before that axis existed plots as what it is.

A `tradeoff` figure used to plot the same runs as points in wall time against
sensitivity. The points overlapped, each carried two parameters, and no reading
of them survived the overlap, so it was deleted.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

# a big canvas with the text scaled up to match, so a figure dropped into a
# document or a slide is legible without anyone zooming
SCALE = 1.75
mpl.rcParams.update({"font.size": mpl.rcParams["font.size"] * SCALE})

PANELS = (22, 9)
SIX = (30, 17)

# small enough to fit inside a heatmap cell, which the base size does not
CELL_SIZE = 11

# the palette the pid and mgnify benchmarks use, so figures from all three sit
# together
TOL_RED = "#CC3311"

HMMER_COLOR = TOL_RED


def save(fig, out, name):
    """Vector, so the figures hold up at whatever size they end up printed at.

    savefig picks its canvas from the extension, so the Agg backend set above
    is only what keeps this from wanting a display.
    """
    path = out / f"{name}.pdf"
    fig.savefig(path)
    plt.close(fig)
    return path


def colorbar(fig, mappable, ax, **kw):
    """A colorbar that stays vector.

    Matplotlib rasterizes the coloured band by default, which is why a figure
    with nothing else raster in it still ships an embedded image. Turning that
    off leaves quads, and giving them a face-coloured edge covers the hairlines
    that otherwise show between them.
    """
    cb = fig.colorbar(mappable, ax=ax, **kw)
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor("face")
    return cb


@dataclass
class Run:
    """One row: a named run, what it found, and what it cost."""

    name: str
    tool: str
    params: dict
    wall_s: float
    hits: int
    sens: float

    def param(self, key):
        """A swept setting as a number, or nan where this run has no such
        setting -- which is what keeps the unpruned run off the grid."""
        try:
            return float(self.params[key])
        except (KeyError, ValueError):
            return float("nan")

    @property
    def a(self):
        return self.param("A")

    @property
    def b(self):
        return self.param("B")

    @property
    def attempts(self):
        return self.param("attempts")


@dataclass
class Table:
    runs: list
    query_count: int
    target_count: int
    hmmer_hits: int

    @property
    def searches(self):
        """Everything but hmmer, which is the thing they are measured
        against rather than one of them."""
        return [r for r in self.runs if r.tool != "hmmer"]

    @property
    def cells(self):
        """The runs that sit on the (A, B) grid."""
        return [r for r in self.searches if np.isfinite(r.a) and np.isfinite(r.b)]

    @property
    def full(self):
        """The unpruned run, which is the ceiling the cells are read against.

        It is a search like any other and sits in the table like one; what
        keeps it off the grid is that `--full-dp` records no -A and no -B.
        """
        off = [c for c in self.searches if not (np.isfinite(c.a) and np.isfinite(c.b))]
        return off[0] if off else None

    @property
    def alphas(self):
        return sorted({c.a for c in self.cells})

    @property
    def betas(self):
        return sorted({c.b for c in self.cells})

    @property
    def attempts(self):
        """The -a the grid was run at, ascending.

        Empty for a run made before -a was swept, which is what keeps those
        summaries drawing as the single surface they are.
        """
        return sorted({c.attempts for c in self.cells if np.isfinite(c.attempts)})

    def surface(self, field, attempts=None, corner=True):
        """The grid as a 2d array, alpha down the rows and beta across.

        One row and one column wider than the grid, holding nothing but the
        unpruned run in the far corner. It records no -A and no -B, so it has
        no cell of its own; the corner is where the surface is heading as the
        pruning relaxes, which is the only place it can be read against.

        `corner` is off for a difference of two surfaces: the unpruned run is
        one run at neither -a, so it cancels to a zero that reads as a measured
        result rather than as an absence.
        """
        alphas, betas = self.alphas, self.betas
        out = np.full((len(alphas) + 1, len(betas) + 1), np.nan)

        for c in self.cells:
            if attempts is not None and c.attempts != attempts:
                continue
            out[alphas.index(c.a), betas.index(c.b)] = getattr(c, field)

        if corner and self.full is not None:
            out[len(alphas), len(betas)] = getattr(self.full, field)

        return out


# what every summary.tbl carries, whatever it swept. everything else in the
# header is a setting, and becomes a param on the run.
FIXED = ("name", "tool", "wall_s", "found", "hits", "sens", "hits_sd", "sens_sd")


def read(path):
    """summary.tbl: `#` lines are metadata and the header, the rest are rows.

    The columns between `tool` and `wall_s` are whatever the pipeline swept, so
    they are read off the header rather than known in advance.
    """
    meta, names, runs = {}, None, []

    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            f = line.lstrip("#").split()
            # `query 200 families 38061 residues ...` and friends
            if len(f) >= 2 and f[0] in ("query", "target", "pairs", "hmmer", "seed"):
                meta[f[0]] = f[1:]
            elif names is None and f and f[0] == "name":
                names = f
            continue

        if not line.split():
            continue

        if names is None:
            raise SystemExit(f"no header in {path}")

        row = dict(zip(names, line.split()))
        runs.append(
            Run(
                name=row["name"],
                tool=row["tool"],
                params={k: v for k, v in row.items() if k not in FIXED and v != "-"},
                wall_s=float(row["wall_s"]),
                hits=int(row["hits"]),
                sens=float(row["sens"]),
            )
        )

    if not runs:
        raise SystemExit(f"no runs in {path}")

    def num(key, i=0):
        try:
            return float(meta[key][i])
        except (KeyError, IndexError, ValueError):
            return float("nan")

    return Table(
        runs=runs,
        query_count=int(num("query")),
        target_count=int(num("target")),
        hmmer_hits=int(num("hmmer")),
    )


def subtitle(g):
    return (
        f"{g.query_count:,} profiles x {g.target_count:,} sequences"
        f"   |   hmmer found {g.hmmer_hits:,}"
    )


# ------------------------------------------------------------------ figures


# a field, the panel title it gets, and how its numbers are written
PANELS_BY_FIELD = [
    ("sens", "sensitivity vs hmmer", "viridis", "{:.2f}"),
    ("wall_s", "wall time (s)", "magma_r", "{:.2f}"),
]

# where nail ships
DEFAULT_CELL = (10.0, 16.0)


def panel(ax, g, data, cmap, fmt, title, norm):
    """One surface: the quads, the ticks, a number in every cell."""
    alphas, betas = g.alphas, g.betas

    # cell edges either side of each integer, so the centres stay on 0..n-1 and
    # the ticks and the default-cell box can be placed by index
    edges_x = np.arange(len(betas) + 2) - 0.5
    edges_y = np.arange(len(alphas) + 2) - 0.5

    # pcolormesh rather than imshow: imshow resamples the grid into a bitmap
    # and embeds that, which is what made the pdf pixelate. this draws one quad
    # per cell and stays vector at any zoom. the edge has to actually be drawn
    # for edgecolors="face" to do its job -- at zero width the quads keep the
    # hairline seams between them.
    im = ax.pcolormesh(
        edges_x, edges_y, data, cmap=cmap, norm=norm,
        edgecolors="face", linewidth=0.4, rasterized=False,
    )
    ax.set_ylim(edges_y[0], edges_y[-1])

    # the last tick on each axis is the unpruned run, which is what -A and -B
    # are approaching rather than a value either of them takes
    ax.set_xticks(range(len(betas) + 1), [f"{b:g}" for b in betas] + ["full"])
    ax.set_yticks(range(len(alphas) + 1), [f"{a:g}" for a in alphas] + ["full"])
    ax.set_xlabel("-B   (global pruning)")
    ax.set_ylabel("-A   (local pruning)")
    ax.set_title(title)

    # a number in every cell: the surface is small enough to read
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            if np.isnan(v):
                continue
            # against the colour actually painted, so the text stays legible
            # whichever end of the map the cell landed on. the weights are
            # relative luminance: green carries most of the brightness the eye
            # sees, and an unweighted mean calls mid magma light when it is not
            r, gr, bl = im.cmap(im.norm(v))[:3]
            light = 0.299 * r + 0.587 * gr + 0.114 * bl
            ax.text(
                j, i, fmt.format(v), ha="center", va="center",
                fontsize=CELL_SIZE, color="black" if light > 0.55 else "white",
            )

    a, b = DEFAULT_CELL
    if a in alphas and b in betas:
        ax.add_patch(
            plt.Rectangle(
                (betas.index(b) - 0.5, alphas.index(a) - 0.5), 1, 1,
                fill=False, edgecolor=HMMER_COLOR, linewidth=3.5,
            )
        )

    return im


def heatmaps(g, out):
    """The surfaces, annotated. Wants an (A, B) grid.

    Two panels for one -a, and six for two: each field at each -a, then the
    difference between them. The difference is the panel the second arm was run
    for -- it says where on the grid nail's disjoint-cloud recovery is doing
    the work, which is sensitivity that -A and -B are credited with and did not
    earn.
    """
    if not g.cells:
        return None

    arms = g.attempts
    if len(arms) < 2:
        return one_arm(g, out)

    lo, hi = arms[0], arms[-1]
    if len(arms) > 2:
        raise SystemExit(f"{len(arms)} -a values in the table; the figure draws two")

    fig, axes = plt.subplots(2, 3, figsize=SIX, constrained_layout=True)

    for row, (field, title, cmap, fmt) in zip(axes, PANELS_BY_FIELD):
        at_hi = g.surface(field, hi)
        at_lo = g.surface(field, lo)
        # the arms share a scale, so the two panels can be read against each
        # other rather than each against itself
        both = Normalize(np.nanmin([at_lo, at_hi]), np.nanmax([at_lo, at_hi]))

        panel(row[0], g, at_hi, cmap, fmt, f"{title}   -a {hi:g}", both)
        im = panel(row[1], g, at_lo, cmap, fmt, f"{title}   -a {lo:g}", both)
        # one bar for the pair, since the point of the shared scale is that
        # there is only one scale to read
        colorbar(fig, im, [row[0], row[1]], shrink=0.9)

        # the corner is dropped rather than differenced: the unpruned run is
        # one run at neither -a, so subtracting it from itself would draw a
        # zero that reads as a measured result
        diff = g.surface(field, hi, corner=False) - g.surface(field, lo, corner=False)
        span = np.nanmax(np.abs(diff)) or 1.0
        im = panel(
            row[2], g, diff, "RdBu_r", fmt,
            f"{title}   -a {hi:g} minus -a {lo:g}",
            Normalize(-span, span),
        )
        colorbar(fig, im, row[2], shrink=0.9)

    fig.suptitle(
        f"cloud search pruning: what -A and -B cost, and what -a recovers\n"
        f"{subtitle(g)}   |   red box is nail's default (-A 10 -B 16)"
    )
    return save(fig, out, "heatmaps")


def one_arm(g, out):
    """Both fields side by side, for a table with a single -a in it."""
    fig, axes = plt.subplots(1, 2, figsize=PANELS, constrained_layout=True)

    for ax, (field, title, cmap, fmt) in zip(axes, PANELS_BY_FIELD):
        data = g.surface(field)
        im = panel(
            ax, g, data, cmap, fmt, title,
            Normalize(np.nanmin(data), np.nanmax(data)),
        )
        colorbar(fig, im, ax, shrink=0.9)

    fig.suptitle(
        f"cloud search pruning: what -A and -B cost\n"
        f"{subtitle(g)}   |   red box is nail's default (-A 10 -B 16)"
    )
    return save(fig, out, "heatmaps")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("summary", help="the summary.tbl to plot")
    ap.add_argument("--out", default="figures", help="where the pdfs go")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # a table with no grid in it draws nothing and says so, rather than failing
    # the run: what a pipeline can be plotted as is a property of what it swept
    g = read(args.summary)
    drawn = [path for path in (heatmaps(g, out),) if path]

    for path in drawn:
        print(f"wrote {path}")

    if not drawn:
        print("skipped: no (A, B) grid to draw")


if __name__ == "__main__":
    main()
