#!/usr/bin/env python3
import argparse
from pathlib import Path

from plot import axes, Scatter, TOL_CYAN, TOL_RED

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, LogNorm


import numpy as np

xmin = None
xmax = None
ymin = None
ymax = None


def heatmap(points, ax):
    N_BINS = 1000
    global xmin, xmax, ymin, ymax

    xbins = np.linspace(xmin, xmax, N_BINS)
    ybins = np.linspace(ymin, ymax, N_BINS)

    def plot(points, cmap):
        colors = cmap(np.linspace(0.2, 1, 256))
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
            norm=LogNorm(),
            shading='auto',
            rasterized=True,
        )

    return plot(points, cm.Blues)


def main(args):
    with open(args.score) as f:
        points = Scatter(f.readlines())

    global xmin, xmax, ymin, ymax
    xmin = min(points.x)
    xmax = max(points.x)
    ymin = min(points.y)
    ymax = max(points.y)

    fig, ax = axes()

    mesh = heatmap(points, ax)

    fig.colorbar(mesh, ax=ax, label="count")
    # ax.legend()

    ax.set_title("Sequence bitscore")
    ax.grid(True)

    ax.set_xlabel("Sequence bitscore from full F/B")
    ax.set_xlim(0, 500)

    ax.set_ylabel("Sequence bitscore from sparse F/B")
    ax.set_ylim(0, 500)

    plt.savefig(args.out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("score", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    main(args)
