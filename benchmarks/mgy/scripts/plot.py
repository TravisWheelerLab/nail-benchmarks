#!/usr/bin/env python3
"""The pruning heatmaps, out of the summary.tbl `mgy parse summary` writes.

    plot.py summary.tbl --out figures/

One row per run, and the settings a run was swept over are columns. This wants
-A and -B among them, so it draws for cloud-search and skips itself for the
pipelines that swept something else.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

# what util/scripts/plot.py does for the other benchmarks: a big canvas with
# the text scaled up to match, so a figure dropped into a document or a slide
# is legible without anyone zooming
SCALE = 1.75
mpl.rcParams.update({"font.size": mpl.rcParams["font.size"] * SCALE})

PANELS = (22, 9)

# small enough to fit inside a heatmap cell, which the base size does not
CELL_SIZE = 11

# the palette pct-id and mgnify use, so figures from all three sit together
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
        seed_wall_s=num("seed", 1),
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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("summary", help="the summary.tbl to plot")
    ap.add_argument("--out", default="figures", help="where the pdfs go")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # a table with no grid in it draws nothing and says so, rather than failing
    # the run: what a pipeline can be plotted as is a property of what it swept
    path = heatmaps(read(args.summary), out)
    print(f"wrote {path}" if path else "skipped: no (A, B) grid to draw")


if __name__ == "__main__":
    main()
