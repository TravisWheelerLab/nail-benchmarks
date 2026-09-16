#!/usr/bin/env python3
"""The cloud-search figures, out of the summary.tbl `mgy parse summary` writes.

    plot.py summary.tbl --out figures/

One row per run, and the settings a run was swept over are columns. These want
-A and -B among them, so they draw for cloud-search and skip themselves for the
pipelines that swept something else.

Two figures. `heatmaps` is the surface -- sensitivity and wall time over the
(A, B) grid. `tradeoff` is the same runs as points in time against
sensitivity, which is where the unpruned `--full-dp` cell earns its keep: it
is the ceiling every pruned cell is measured against, and it sits off the grid
because it has no -A or -B of its own.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.ticker import NullFormatter, NullLocator

# a big canvas with the text scaled up to match, so a figure dropped into a
# document or a slide is legible without anyone zooming
SCALE = 1.75
mpl.rcParams.update({"font.size": mpl.rcParams["font.size"] * SCALE})

WIDE = (16, 9)
PANELS = (22, 9)

# small enough to fit inside a heatmap cell, which the base size does not
CELL_SIZE = 11

# the palette the pid and mgnify benchmarks use, so figures from all three sit
# together
TOL_RED = "#CC3311"
TOL_TEAL = "#009988"

HMMER_COLOR = TOL_RED
FULL_COLOR = TOL_TEAL


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


def log_time(ax):
    """A log time axis that keeps its labels.

    The grid spans well under a decade, so the default locator finds one or two
    ticks and the axis reads as blank. Asking for the 1/2/3/5/7 subdivisions
    fills it in.
    """
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(plt.LogLocator(base=10, subs=(1, 2, 3, 5, 7)))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())


def pareto(points):
    """The points nothing else beats on both time and sensitivity."""
    front, best = [], -np.inf
    for p in sorted(points, key=lambda p: (p.wall_s, -p.sens)):
        if p.sens > best:
            front.append(p)
            best = p.sens
    return front


def alpha_colors(g):
    cmap = plt.get_cmap("viridis")
    norm = Normalize(min(g.alphas), max(g.alphas))
    return cmap, norm


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


@dataclass
class Table:
    runs: list
    query_count: int
    target_count: int
    hmmer_hits: int
    hmmer_wall_s: float
    seed_wall_s: float

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

    def surface(self, field):
        """The grid as a 2d array, alpha down the rows and beta across."""
        alphas, betas = self.alphas, self.betas
        out = np.full((len(alphas), len(betas)), np.nan)
        for c in self.cells:
            out[alphas.index(c.a), betas.index(c.b)] = getattr(c, field)
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
        hmmer_wall_s=num("hmmer", 2),
        seed_wall_s=num("seed", 0),
    )


def subtitle(g):
    return (
        f"{g.query_count:,} profiles x {g.target_count:,} sequences"
        f"   |   hmmer found {g.hmmer_hits:,}"
    )


# ------------------------------------------------------------------ figures


def heatmaps(g, out):
    """Both surfaces side by side, annotated. Wants an (A, B) grid."""
    if not g.cells:
        return None

    alphas, betas = g.alphas, g.betas
    fig, axes = plt.subplots(1, 2, figsize=PANELS, constrained_layout=True)

    panels = [
        ("sens", "sensitivity vs hmmer", "viridis", "{:.2f}"),
        ("wall_s", "wall time (s)", "magma_r", "{:.2f}"),
    ]

    # cell edges either side of each integer, so the centres stay on 0..n-1 and
    # the ticks and the default-cell box can be placed by index
    edges_x = np.arange(len(betas) + 1) - 0.5
    edges_y = np.arange(len(alphas) + 1) - 0.5

    for ax, (field, title, cmap, fmt) in zip(axes, panels):
        data = g.surface(field)

        # pcolormesh rather than imshow: imshow resamples the grid into a
        # bitmap and embeds that, which is what made the pdf pixelate. this
        # draws one quad per cell and stays vector at any zoom. the edge has to
        # actually be drawn for edgecolors="face" to do its job -- at zero
        # width the quads keep the hairline seams between them.
        im = ax.pcolormesh(
            edges_x, edges_y, data, cmap=cmap,
            edgecolors="face", linewidth=0.4, rasterized=False,
        )
        ax.set_ylim(edges_y[0], edges_y[-1])

        ax.set_xticks(range(len(betas)), [f"{b:g}" for b in betas])
        ax.set_yticks(range(len(alphas)), [f"{a:g}" for a in alphas])
        ax.set_xlabel("-B   (global pruning)")
        ax.set_ylabel("-A   (local pruning)")
        ax.set_title(title)

        # a number in every cell: the surface is small enough to read
        norm = Normalize(np.nanmin(data), np.nanmax(data))
        for i in range(len(alphas)):
            for j in range(len(betas)):
                v = data[i, j]
                if np.isnan(v):
                    continue
                shade = "white" if norm(v) < 0.55 else "black"
                if cmap.endswith("_r"):
                    shade = "black" if norm(v) < 0.45 else "white"
                ax.text(
                    j, i, fmt.format(v), ha="center", va="center",
                    fontsize=CELL_SIZE, color=shade,
                )

        colorbar(fig, im, ax, shrink=0.9)

        # where nail ships
        if 10.0 in alphas and 16.0 in betas:
            ax.add_patch(
                plt.Rectangle(
                    (betas.index(16.0) - 0.5, alphas.index(10.0) - 0.5), 1, 1,
                    fill=False, edgecolor=HMMER_COLOR, linewidth=3.5,
                )
            )

    fig.suptitle(
        f"cloud search pruning: what -A and -B cost\n"
        f"{subtitle(g)}   |   red box is nail's default (-A 10 -B 16)"
    )
    return save(fig, out, "heatmaps")


def tradeoff(g, out):
    """Every cell as a point in time/sensitivity, with the front traced."""
    if not g.cells:
        return None

    fig, ax = plt.subplots(figsize=WIDE, constrained_layout=True)
    cmap, norm = alpha_colors(g)

    for c in g.cells:
        ax.scatter(
            c.wall_s, c.sens, s=190, color=cmap(norm(c.a)),
            edgecolor="white", linewidth=1.0, zorder=3,
        )

    front = pareto(g.cells)
    ax.plot(
        [p.wall_s for p in front], [p.sens for p in front],
        color="black", linewidth=2.0, linestyle="--", zorder=2,
        label="pareto front",
    )

    # the ceiling: the most nail can find off these seeds, and the longest it
    # can take to find it
    if g.full is not None:
        ax.axhline(
            g.full.sens, color=FULL_COLOR, linewidth=2.5, linestyle=":",
            label=f"--full-dp ceiling ({g.full.sens:.3f})",
        )
        ax.scatter(
            [g.full.wall_s], [g.full.sens], marker="*", s=900,
            color=FULL_COLOR, edgecolor="white", linewidth=1.2, zorder=4,
            label=f"--full-dp ({g.full.wall_s:.2f}s)",
        )

    if np.isfinite(g.hmmer_wall_s):
        ax.axvline(
            g.hmmer_wall_s, color=HMMER_COLOR, linewidth=2.5, linestyle="-.",
            label=f"hmmer wall ({g.hmmer_wall_s:.2f}s)",
        )

    default = next((c for c in g.cells if c.a == 10.0 and c.b == 16.0), None)
    if default is not None:
        ax.annotate(
            "nail default\n-A 10 -B 16",
            (default.wall_s, default.sens),
            textcoords="offset points", xytext=(30, -70), fontweight="bold",
            arrowprops=dict(arrowstyle="->", color="black", linewidth=1.8),
        )

    colorbar(fig, plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax,
             label="-A   (local pruning)")

    # most of the grid sits in the cheap corner, so a linear axis stacks it all
    # against the left edge
    log_time(ax)

    ax.set_xlabel("wall time (s), log scale")
    ax.set_ylabel("sensitivity vs hmmer")
    ax.set_title(f"what pruning buys and what it costs\n{subtitle(g)}")
    ax.grid(alpha=0.25, zorder=0, which="both")
    ax.legend(loc="lower right")

    return save(fig, out, "tradeoff")


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
    drawn = [path for path in (heatmaps(g, out), tradeoff(g, out)) if path]

    for path in drawn:
        print(f"wrote {path}")

    if not drawn:
        print("skipped: no (A, B) grid to draw")


if __name__ == "__main__":
    main()
