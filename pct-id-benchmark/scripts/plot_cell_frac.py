#!/usr/bin/env python3
import argparse

from plot import Scatter, TOL_CYAN, TOL_MAGENTA

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap


import numpy as np

xmin = None
xmax = None
ymin = None
ymax = None


def heatmap(true_points, decoy_points, ax, long_points=None):
    N_BINS = 200
    global xmin, xmax, ymin, ymax

    xbins = np.logspace(np.log10(xmin), np.log10(xmax), N_BINS)
    ybins = np.logspace(np.log10(ymin), np.log10(ymax), N_BINS)

    def plot(points, cmap):
        colors = cmap(np.linspace(0.25, 1, 256))
        colors[:, -1] = np.linspace(0.8, 1, 256)
        cmap = ListedColormap(colors)
        h, xedges, yedges = np.histogram2d(
            points.x,
            points.y,
            bins=[xbins, ybins]
        )

        h[h == 0] = np.nan

        ax.pcolormesh(
            xedges,
            yedges,
            h.T,
            cmap=cmap,
            shading='auto'
        )

    plot(decoy_points, cm.Reds)
    plot(true_points, cm.Blues)

    if long_points:
        ax.scatter(
            long_points.x, long_points.y,
            label="Long sequences",
            color=TOL_CYAN,
            marker='*',
        )


def scatter(true_points, decoy_points, ax):
    alpha = 0.6

    ax.scatter(
        decoy_points.x, decoy_points.y,
        label="decoys",
        color=TOL_MAGENTA,
        marker='v',
        alpha=alpha,
        edgecolors="none",
        rasterized=True,
    )

    ax.scatter(
        true_points.x, true_points.y,
        label="True positives",
        color=TOL_CYAN,
        marker='^',
        alpha=alpha,
        edgecolors="none",
        rasterized=True,
    )

    ax.legend()


def main(args):
    with open(args.true_hits) as f:
        true_points = Scatter(f.readlines())

    with open(args.decoy_hits) as f:
        decoy_points = Scatter(f.readlines())

    global xmin, xmax, ymin, ymax
    xmin = min(min(true_points.x), min(decoy_points.x))
    xmax = max(max(true_points.x), max(decoy_points.x))
    ymin = min(min(true_points.y), min(decoy_points.y))
    ymax = max(max(true_points.y), max(decoy_points.y))

    long_points = None
    if args.long_hits:
        with open(args.long_hits) as f:
            long_points = Scatter(f.readlines())

    fig, ax = plt.subplots(figsize=(16, 9))

    # scatter(true_points, decoy_points, ax)
    heatmap(true_points, decoy_points, ax, long_points)

    ax.set_title("Cells computed")
    ax.grid(True)

    ax.set_xlabel("Total DP cells (product of target and query lengths)")
    ax.set_xscale("log")
    ax.set_xlim(1e4, 1e10)

    ax.set_ylabel("Fraction of DP cells computed by sparse F/B")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1)

    plt.savefig(args.true_hits.with_name("cells.pdf"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("true_hits", type=Path)
    p.add_argument("decoy_hits", type=Path)
    p.add_argument("--long_hits", type=Path)
    args = p.parse_args()
    main(args)
