#!/usr/bin/env python3
import argparse

from plot import axes, Scatter, TOL_CYAN, TOL_RED, TOL_ORANGE

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap


import numpy as np

xmin = None
xmax = None
ymin = None
ymax = None


def heatmap(points, ax):
    N_BINS = 500
    global xmin, xmax, ymin, ymax

    xbins = np.linspace(xmin, xmax, N_BINS)
    ybins = np.linspace(ymin, ymax, N_BINS)

    def plot(points, cmap):
        colors = cmap(np.linspace(0.5, 1, 256))
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

    plot(points, cm.Blues)


def scatter(points, ax):
    alpha = 0.6

    ax.scatter(
        points.x, points.y,
        label="true positives",
        color=TOL_CYAN,
        # marker='^',
        alpha=alpha,
        edgecolors="none",
        rasterized=True,
    )


def yx(ax):
    # y = x
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    lo = min(xmin, ymin)
    hi = max(xmax, ymax)
    ax.plot(
        [lo, hi], [lo, hi],
        color=TOL_ORANGE,
        linestyle="--",
        linewidth=1,
        label="y=x"
    )


def trend(points, ax):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    lo = min(xmin, ymin)
    hi = max(xmax, ymax)

    m, b = np.polyfit(points.x, points.y, 1)
    x_fit = np.array([lo, hi])
    y_fit = m * x_fit + b
    ax.plot(
        x_fit, y_fit,
        color=TOL_RED,
        linewidth=1,
        label="trend"
    )


def main(args):
    with open(args.score) as f:
        points = Scatter(f.readlines())

    global xmin, xmax, ymin, ymax
    xmin = min(points.x)
    xmax = max(points.x)
    ymin = min(points.y)
    ymax = max(points.y)

    fig, ax = axes()

    scatter(points, ax)
    # heatmap(points, ax)
    yx(ax)
    trend(points, ax)

    ax.legend()

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
