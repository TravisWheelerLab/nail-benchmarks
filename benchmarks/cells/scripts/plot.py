#!/usr/bin/env python3
import argparse
from pathlib import Path

from collections.abc import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.legend_handler import HandlerBase

import numpy as np

TOL_CYAN = "#33BBEE"

# drawn in this order, so the original is on top
FORMS = [
    ("shuffled", "Shuffled", cm.Greens),
    ("reversed", "Reversed", cm.Reds),
    ("original", "Original", cm.Blues),
]

xmin = None
xmax = None
ymin = None
ymax = None


def axes():
    scale = 1.75
    mpl.rcParams.update({
        "font.size": mpl.rcParams["font.size"] * scale,
    })

    fig, ax = plt.subplots(
        figsize=(16, 9),
        # constrained_layout=True
    )

    P = 0.075
    fig.subplots_adjust(
        left=P,
        right=1.0 - P,
        bottom=P,
        top=1.0 - P
    )

    return fig, ax


class Scatter:
    x: [float]
    y: [float]

    def __init__(self, lines: Iterable[str]) -> None:
        """The first two columns of a `#`-headed table, as x and y."""
        self.x = []
        self.y = []
        for line in lines:
            if not line.strip() or line.startswith("#"):
                continue
            x, y = line.split()[:2]
            self.x.append(float(x))
            self.y.append(float(y))


class HandlerColormap(HandlerBase):
    def __init__(self, cmap, norm, n=256):
        self.cmap = cmap
        self.norm = norm
        self.n = n
        super().__init__()

    def create_artists(self, legend, orig_handle,
                       x0, y0, w, h, fontsize, trans):
        artists = []

        # gradient fills full handle area so it centers with the label text
        for i in range(self.n):
            xi = x0 + w * i / self.n
            wi = w / self.n
            artists.append(mpl.patches.Rectangle(
                (xi, y0), wi, h,
                transform=trans,
                facecolor=self.cmap(i / self.n),
                edgecolor='none',
            ))

        # tick marks and labels below the gradient (in labelspacing gap)
        vmin = self.norm.vmin
        vmax = self.norm.vmax
        if vmin is not None and vmax is not None:
            for val in np.linspace(vmin, vmax, 5):
                frac = (val - vmin) / (vmax - vmin)
                xt = x0 + frac * w
                tick_bot = y0 - fontsize * 0.15
                artists.append(mpl.lines.Line2D(
                    [xt, xt], [y0, tick_bot],
                    transform=trans,
                    color='black',
                    linewidth=0.8,
                ))
                artists.append(mpl.text.Text(
                    xt, y0 - fontsize * 0.4,
                    str(int(val)),
                    transform=trans,
                    fontsize=fontsize * 0.55,
                    ha='center', va='top',
                    clip_on=False,
                ))

        return artists


def heatmap(forms, ax, long_points=None):
    N_BINS = 200
    global xmin, xmax, ymin, ymax

    xbins = np.logspace(np.log10(xmin), np.log10(xmax), N_BINS)
    ybins = np.logspace(np.log10(ymin), np.log10(ymax), N_BINS)

    def plot(points, cmap_base):
        colors = cmap_base(np.linspace(0.25, 1, 256))
        colors[:, -1] = np.linspace(0.8, 1, 256)
        cmap = ListedColormap(colors)
        h, xedges, yedges = np.histogram2d(
            points.x,
            points.y,
            bins=[xbins, ybins]
        )
        h[h == 0] = np.nan
        return ax.pcolormesh(
            xedges,
            yedges,
            h.T,
            cmap=cmap,
            shading='auto',
            rasterized=True,
        )

    handles = []
    labels = []
    handler_map = {}

    for points, label, cmap_base in forms:
        mesh = plot(points, cmap_base)
        mesh.set_norm(mpl.colors.Normalize(
            vmin=np.nanmin(mesh.get_array()),
            vmax=np.nanmax(mesh.get_array()),
        ))

        sm = plt.cm.ScalarMappable(norm=mesh.norm, cmap=mesh.cmap)
        handles.append(sm)
        labels.append(label)
        handler_map[sm] = HandlerColormap(mesh.cmap, mesh.norm)

    if long_points:
        handles.append(ax.scatter(
            long_points.x, long_points.y,
            color=TOL_CYAN,
            marker='^',
            s=64,
        ))
        labels.append("Long sequences")

    ax.legend(
        handles=handles,
        labels=labels,
        handler_map=handler_map,
        handleheight=1.0,
        handlelength=6.0,
        labelspacing=1.2,
    )


def main(args):
    forms = []
    for form, label, cmap_base in FORMS:
        with open(getattr(args, form)) as f:
            points = Scatter(f.readlines())
        if not points.x:
            print(f"warning: no {form} points to draw")
            continue
        zeros = sum(1 for y in points.y if y <= 0)
        if zeros:
            print(f"warning: {zeros} {form} points have a fraction of 0, "
                  "which a log axis cannot draw")
        forms.append((points, label, cmap_base))

    global xmin, xmax, ymin, ymax
    xmin = min(min(p.x) for p, _, _ in forms)
    xmax = max(max(p.x) for p, _, _ in forms)
    ymin = min(min(y for y in p.y if y > 0) for p, _, _ in forms)
    ymax = max(max(p.y) for p, _, _ in forms)

    long_points = None
    if args.long_hits:
        with open(args.long_hits) as f:
            long_points = Scatter(f.readlines())

    fig, ax = axes()

    heatmap(forms, ax, long_points)

    ax.set_title("Cells computed")
    ax.grid(True)

    ax.set_xlabel("Total DP cells (product of target and query lengths)")
    ax.set_xscale("log")
    ax.set_xlim(1e4, 1e10)

    ax.set_ylabel("Fraction of DP cells computed by sparse F/B")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1)

    plt.savefig(args.out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for form, _, _ in FORMS:
        p.add_argument(f"--{form}", type=Path, required=True)
    p.add_argument("out", type=Path)
    p.add_argument("--long_hits", type=Path)
    args = p.parse_args()
    main(args)
