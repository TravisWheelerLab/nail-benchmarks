#!/usr/bin/env python3
"""The cloud-search figures, out of the summary.tbl `cloud-search parse summary` writes.

    plot.py summary.tbl --out figures/

One row per run, and the settings a run was swept over are columns. These want
-A and -B among them, so they draw for cloud-search and skip themselves for the
pipelines that swept something else.

One figure. `cloud-search` is the surface over the (A, B) grid, with the unpruned
`--full-dp` run in the far corner. A cell's colour is its coverage, its hits
over nail's cutoff as a percentage of the full run's, and its text is its
median wall time over the reps. The full run records no -A and no -B, so it
has no cell of its own; the corner is where the surface is heading as the
pruning relaxes, and it is the ceiling every pruned cell is read against.

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
from matplotlib.colors import LinearSegmentedColormap, LogNorm

# a big canvas with the text scaled up to match, so a figure dropped into a
# document or a slide is legible without anyone zooming
SCALE = 1.75
mpl.rcParams.update({"font.size": mpl.rcParams["font.size"] * SCALE})

# pid's canvas, so the figures sit together
FIGSIZE = (16, 9)

# small enough to fit inside a heatmap cell, which the base size does not
CELL_SIZE = 13

# pid's palette, where nail is magenta. coverage runs from pale through
# nail's magenta to a dark purple, several hues so the top few percent still
# read apart, and the default cell is boxed in the palette's blue
TOL_MAGENTA = "#EE3377"
TOL_BLUE = "#0077BB"
NAIL_MAP = LinearSegmentedColormap.from_list(
    "nail", ["#FDF0E6", "#F6B48F", TOL_MAGENTA, "#7B1F6A", "#2B0A3D"]
)

DEFAULT_COLOR = TOL_BLUE


def save(fig, out, name):
    """Vector, so the figures hold up at whatever size they end up printed at.

    savefig picks its canvas from the extension, so the Agg backend set above
    is only what keeps this from wanting a display.
    """
    path = out / f"{name}.pdf"
    # cropped to what is drawn: the grid is square, so a fixed canvas leaves
    # blank bands either side of it
    fig.savefig(path, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    return path


def colorbar(fig, mappable, ax, **kw):
    """A colorbar that stays vector.

    Matplotlib rasterizes the coloured band by default, which is why a figure
    with nothing else raster in it still ships an embedded image. Turning that
    off leaves quads, and giving them a face-coloured edge covers the hairlines
    that otherwise show between them.
    """
    cb = fig.colorbar(mappable, ax=ax, **kw) if ax is not None else fig.colorbar(mappable, **kw)
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor("face")
    return cb


@dataclass
class Run:
    """One row: a named run, what it found, and what it cost."""

    name: str
    params: dict
    reps: int
    wall_s: float
    hits: int
    coverage: float

    @property
    def percent(self):
        return 100 * self.coverage

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


@dataclass
class Table:
    runs: list
    query_count: int
    target_count: int
    full_hits: int

    @property
    def cells(self):
        """The runs that sit on the (A, B) grid."""
        return [r for r in self.runs if np.isfinite(r.a) and np.isfinite(r.b)]

    @property
    def full(self):
        """The unpruned run, which is the ceiling the cells are read against.

        It is a search like any other and sits in the table like one; what
        keeps it off the grid is that `--full-dp` records no -A and no -B.
        """
        off = [c for c in self.runs if not (np.isfinite(c.a) and np.isfinite(c.b))]
        return off[0] if off else None

    @property
    def alphas(self):
        return sorted({c.a for c in self.cells})

    @property
    def betas(self):
        return sorted({c.b for c in self.cells})

    def surface(self, field):
        """The grid as a 2d array, alpha down the rows and beta across.

        One row and one column wider than the grid, holding nothing but the
        unpruned run in the far corner. It records no -A and no -B, so it has
        no cell of its own; the corner is where the surface is heading as the
        pruning relaxes, which is the only place it can be read against.
        """
        alphas, betas = self.alphas, self.betas
        out = np.full((len(alphas) + 1, len(betas) + 1), np.nan)

        for c in self.cells:
            out[alphas.index(c.a), betas.index(c.b)] = getattr(c, field)

        if self.full is not None:
            out[len(alphas), len(betas)] = getattr(self.full, field)

        return out


# what every summary.tbl carries, whatever it swept. everything else in the
# header is a setting, and becomes a param on the run.
FIXED = (
    "name",
    "reps",
    "wall_s",
    "wall_min",
    "wall_max",
    "hits",
    "coverage",
    "left_out",
)


def read(path):
    """summary.tbl: `#` lines are metadata and the header, the rest are rows.

    The columns between `name` and `wall_s` are whatever the pipeline swept, so
    they are read off the header rather than known in advance.
    """
    meta, names, runs = {}, None, []

    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            f = line.lstrip("#").split()
            # `query 200 families 38061 residues ...` and friends
            if len(f) >= 2 and f[0] in ("query", "target", "pairs", "full", "seed"):
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
                params={k: v for k, v in row.items() if k not in FIXED and v != "-"},
                reps=int(row["reps"]),
                wall_s=float(row["wall_s"]),
                hits=int(row["hits"]),
                coverage=float(row["coverage"]),
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
        full_hits=int(num("full")),
    )


# ------------------------------------------------------------------ figures


# where nail ships
DEFAULT_CELL = (10.0, 16.0)


def wall(v):
    """A wall time short enough for a cell: three significant figures or so."""
    return (f"{v:.2f}" if v < 10 else f"{v:.1f}" if v < 100 else f"{v:.0f}") + "s"


def heatmap(g, out):
    """Coverage as colour and median wall time as text over the (A, B) grid.
    Wants a grid."""
    if not g.cells:
        return None

    alphas, betas = g.alphas, g.betas
    # coloured by the share of the full run's hits a cell misses, on a log
    # scale, so 99.5% and 99.9% read apart as clearly as 60% and 90% do. the
    # full run misses none, so it takes the scale's floor
    colour = 100.0 - g.surface("percent")
    floor = 0.05
    colour = np.where(np.isnan(colour), np.nan, np.maximum(colour, floor))
    text = g.surface("wall_s")

    fig, ax = plt.subplots(figsize=FIGSIZE)
    # the grid is square, so it is centred on the canvas with its colorbar
    # beside it rather than laid out by an engine that pushes it to one side
    fig.subplots_adjust(left=0.075, right=0.925, bottom=0.14, top=0.92)

    # cell edges either side of each integer, so the centres stay on 0..n-1 and
    # the ticks and the default-cell box can be placed by index
    edges_x = np.arange(len(betas) + 2) - 0.5
    edges_y = np.arange(len(alphas) + 2) - 0.5

    # pcolormesh rather than imshow: imshow resamples the grid into a bitmap
    # and embeds that, which is what made the pdf pixelate. this draws one quad
    # per cell and stays vector at any zoom. a white edge separates the cells
    # and covers the hairline seams between quads
    im = ax.pcolormesh(
        edges_x, edges_y, colour, cmap=NAIL_MAP.reversed(),
        norm=LogNorm(floor, max(np.nanmax(colour), 1.0)),
        edgecolors="white", linewidth=1.5, rasterized=False,
    )
    ax.set_ylim(edges_y[0], edges_y[-1])
    ax.set_aspect("equal")
    # a frame and ticks like pid's axes, round the grid and its full row
    for side in ax.spines.values():
        side.set_visible(True)
        side.set_color("black")

    # the last tick on each axis is the unpruned run, which is what -A and -B
    # are approaching rather than a value either of them takes
    ax.set_xticks(range(len(betas) + 1), [f"{b:g}" for b in betas] + ["full"])
    ax.set_yticks(range(len(alphas) + 1), [f"{a:g}" for a in alphas] + ["full"])
    ax.set_xlabel("Global pruning threshold (-B)")
    ax.set_ylabel("Local pruning threshold (-A)")

    for i in range(colour.shape[0]):
        for j in range(colour.shape[1]):
            if np.isnan(colour[i, j]):
                # the full row and column hold one run, in the corner; the
                # rest of them are settings no run took
                if i == len(alphas) or j == len(betas):
                    ax.text(j, i, "\u2013", ha="center", va="center",
                            fontsize=CELL_SIZE * 2, color="0.4")
                continue
            # against the colour actually painted, so the text stays legible
            # whichever end of the map the cell landed on. the weights are
            # relative luminance: green carries most of the brightness the eye
            # sees
            r, gr, bl = im.cmap(im.norm(colour[i, j]))[:3]
            light = 0.299 * r + 0.587 * gr + 0.114 * bl
            ax.text(
                j, i, wall(text[i, j]), ha="center", va="center",
                fontsize=CELL_SIZE, color="black" if light > 0.55 else "white",
            )

    a, b = DEFAULT_CELL
    if a in alphas and b in betas:
        ax.add_patch(
            plt.Rectangle(
                (betas.index(b) - 0.5, alphas.index(a) - 0.5), 1, 1,
                fill=False, edgecolor=DEFAULT_COLOR, linewidth=3.5,
            )
        )
        ax.legend(
            [plt.Rectangle((0, 0), 1, 1, fill=False, edgecolor=DEFAULT_COLOR, linewidth=3)],
            [f"nail default (-A {a:g} -B {b:g})"],
            loc="upper center", bbox_to_anchor=(0.5, -0.1), frameon=False,
            handlelength=1.2, handleheight=1.2,
        )

    cb = colorbar(fig, im, None, cax=ax.inset_axes([1.06, 0.0, 0.035, 0.88]))
    cb.set_label("Coverage (% of hits recovered with --full-dp)")
    # ticks placed by what a cell misses and labelled by what it covers
    misses = [m for m in (0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50) if m <= np.nanmax(colour)]
    cb.set_ticks(misses, labels=[f"{100 - m:g}%" for m in misses])
    cb.ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
    cb.ax.invert_yaxis()
    cb.outline.set_visible(False)

    ax.set_title("Cloud Search: tradeoff between sensitivity and runtime", pad=14)
    return save(fig, out, "cloud-search")


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
    drawn = [path for path in (heatmap(g, out),) if path]

    for path in drawn:
        print(f"wrote {path}")

    if not drawn:
        print("skipped: no (A, B) grid to draw")


if __name__ == "__main__":
    main()
