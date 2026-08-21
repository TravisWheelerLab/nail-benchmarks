#!/usr/bin/env python3
"""Figures for the cloud-search benchmark.

Reads the grid.tbl that `cloud-search parse grid` writes and draws every
figure into an output directory, as pdf.

    plot.py grid.tbl --out figures/
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

# what util/scripts/plot.py does for the other benchmarks: a big canvas with
# the text scaled up to match, so a figure dropped into a document or a slide
# is legible without anyone zooming
SCALE = 1.75
mpl.rcParams.update({"font.size": mpl.rcParams["font.size"] * SCALE})

WIDE = (16, 9)
PANELS = (22, 9)

# small enough to fit inside a heatmap cell, which the base size does not
CELL_SIZE = 11

# the palette pct-id and mgnify use, so figures from all three sit together
TOL_ORANGE = "#EE7733"
TOL_BLUE = "#0077BB"
TOL_CYAN = "#33BBEE"
TOL_TEAL = "#009988"
TOL_RED = "#CC3311"

TOL_VIBRANT = [TOL_ORANGE, TOL_BLUE, TOL_CYAN, TOL_TEAL, TOL_RED]

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


def log_axis(ax, values):
    """Put the x axis on a log scale labelled at the values we actually swept.

    Left alone, matplotlib keeps its own decade and minor labels, which land on
    top of these and turn the axis into mush.
    """
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(plt.FixedLocator(values))
    ax.xaxis.set_major_formatter(plt.FixedFormatter([f"{v:g}" for v in values]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())


@dataclass
class Cell:
    a: float
    b: float
    wall_s: float
    hits: int
    sens: float


@dataclass
class Grid:
    cells: list
    full: Cell
    query_count: int
    target_count: int
    hmmer_hits: int
    hmmer_wall_s: float
    seed_wall_s: float

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


def read(path):
    """grid.tbl: `#` lines are header, the rest are one cell each."""
    cells, full = [], None
    meta = {}

    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            f = line.lstrip("#").split()
            # `query 200 families 38061 residues ...` and friends
            if len(f) >= 2 and f[0] in ("query", "target", "pairs", "hmmer", "seed"):
                meta[f[0]] = f[1:]
            continue

        f = line.split()
        if len(f) != 5:
            continue

        wall, hits, sens = float(f[2]), int(f[3]), float(f[4])
        if f[0] == "-":
            full = Cell(np.nan, np.nan, wall, hits, sens)
        else:
            cells.append(Cell(float(f[0]), float(f[1]), wall, hits, sens))

    if not cells:
        raise SystemExit(f"no cells in {path}")

    def num(key, i=0):
        try:
            return float(meta[key][i])
        except (KeyError, IndexError, ValueError):
            return float("nan")

    return Grid(
        cells=cells,
        full=full,
        query_count=int(num("query")),
        target_count=int(num("target")),
        hmmer_hits=int(num("hmmer")),
        hmmer_wall_s=num("hmmer", 2),
        seed_wall_s=num("seed", 1),
    )


def subtitle(g):
    return (
        f"{g.query_count:,} profiles x {g.target_count:,} sequences"
        f"   |   hmmer found {g.hmmer_hits:,}"
    )


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


# ------------------------------------------------------------------ figures


def heatmaps(g, out):
    """Both surfaces side by side, annotated."""
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
    save(fig, out, "heatmaps")


def tradeoff(g, out):
    """Every cell as a point in time/sensitivity, with the front traced."""
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

    save(fig, out, "tradeoff")


def curves(g, out):
    """One line per -A, so each axis can be read on its own."""
    fig, axes = plt.subplots(1, 2, figsize=PANELS, constrained_layout=True)
    cmap, norm = alpha_colors(g)

    for a in g.alphas:
        row = sorted((c for c in g.cells if c.a == a), key=lambda c: c.b)
        color = cmap(norm(a))
        axes[0].plot(
            [c.b for c in row], [c.sens for c in row],
            marker="o", markersize=9, linewidth=2.2, color=color, label=f"A={a:g}",
        )
        axes[1].plot(
            [c.b for c in row], [c.wall_s for c in row],
            marker="o", markersize=9, linewidth=2.2, color=color,
        )

    for ax, ylab in zip(axes, ["sensitivity vs hmmer", "wall time (s)"]):
        ax.set_xlabel("-B   (global pruning)")
        ax.set_ylabel(ylab)
        log_axis(ax, g.betas)
        ax.grid(alpha=0.25)

    if g.full is not None:
        axes[0].axhline(
            g.full.sens, color=FULL_COLOR, linestyle=":", linewidth=2.5,
            label="--full-dp",
        )
        axes[1].axhline(g.full.wall_s, color=FULL_COLOR, linestyle=":", linewidth=2.5)

    if np.isfinite(g.hmmer_wall_s):
        axes[1].axhline(
            g.hmmer_wall_s, color=HMMER_COLOR, linestyle="-.", linewidth=2.5,
            label="hmmer",
        )
        axes[1].legend(loc="upper left")

    axes[0].legend(ncol=2, loc="lower right", fontsize=CELL_SIZE + 2)
    fig.suptitle(f"each -A across the -B ladder\n{subtitle(g)}")

    save(fig, out, "curves")


def relative(g, out):
    """Everything as a fraction of the unpruned cell, which is the ceiling."""
    if g.full is None:
        return

    fig, ax = plt.subplots(figsize=WIDE, constrained_layout=True)
    cmap, norm = alpha_colors(g)

    for a in g.alphas:
        row = sorted((c for c in g.cells if c.a == a), key=lambda c: c.b)
        ax.plot(
            [c.wall_s / g.full.wall_s for c in row],
            [c.sens / g.full.sens for c in row],
            marker="o", markersize=8, linewidth=2.0,
            color=cmap(norm(a)), label=f"A={a:g}",
        )

    ax.axhline(1.0, color=FULL_COLOR, linestyle=":", linewidth=2.5)
    ax.axvline(1.0, color=FULL_COLOR, linestyle=":", linewidth=2.5)
    ax.scatter([1.0], [1.0], marker="*", s=900, color=FULL_COLOR,
               edgecolor="white", linewidth=1.2, zorder=5, label="--full-dp")

    # the cheapest cell that keeps most of the ceiling, at a few tolerances.
    # one number would be arbitrary; three show how the price climbs as the
    # last percent gets bought back
    for keep, color in zip((0.95, 0.98, 0.99), TOL_VIBRANT[:3]):
        knee = min(
            (c for c in g.cells if c.sens >= keep * g.full.sens),
            key=lambda c: c.wall_s, default=None,
        )
        if knee is None:
            continue
        x = knee.wall_s / g.full.wall_s
        ax.axvline(x, color=color, linewidth=2.2, alpha=0.9, zorder=1)
        ax.annotate(
            f"{keep:.0%} of the ceiling for {x:.0%} of the time"
            f"   (-A {knee.a:g} -B {knee.b:g})",
            (x, 0.40), rotation=90, color=color, fontweight="bold",
            fontsize=CELL_SIZE + 3, ha="right", va="bottom",
        )

    log_time(ax)

    ax.set_xlabel("wall time, as a fraction of --full-dp (log scale)")
    ax.set_ylabel("sensitivity, as a fraction of --full-dp")
    ax.set_title(f"pruning against the unpruned ceiling\n{subtitle(g)}")
    ax.grid(alpha=0.25, which="both")
    ax.legend(ncol=2, loc="lower right", fontsize=CELL_SIZE + 2)

    save(fig, out, "relative")


FIGURES = {
    "heatmaps": heatmaps,
    "tradeoff": tradeoff,
    "curves": curves,
    "relative": relative,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("grid", help="the grid.tbl to plot")
    ap.add_argument("--out", default="figures", help="where the pdfs go")
    ap.add_argument(
        "--only", choices=sorted(FIGURES), action="append",
        help="draw just this figure; repeatable",
    )
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    g = read(args.grid)

    for name in args.only or sorted(FIGURES):
        FIGURES[name](g, out)
        print(f"wrote {out / (name + '.pdf')}")


if __name__ == "__main__":
    main()
