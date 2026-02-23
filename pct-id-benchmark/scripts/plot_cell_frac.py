#!/usr/bin/env python3
import argparse

from plot import Scatter, TOL_CYAN, TOL_MAGENTA

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.legend_handler import HandlerBase


import numpy as np

xmin = None
xmax = None
ymin = None
ymax = None


class HandlerColormap(HandlerBase):
    def __init__(self, cmap, norm, n=256):
        self.cmap = cmap
        self.norm = norm
        self.n = n
        super().__init__()

    def create_artists(self, legend, orig_handle,
                       x0, y0, w, h, fontsize, trans):
        artists = []
        for i in range(self.n):
            xi = x0 + w * i / self.n
            wi = w / self.n
            color = self.cmap(i / self.n)
            r = mpl.patches.Rectangle(
                (xi, y0), wi, h,
                transform=trans,
                facecolor=color,
                edgecolor='none'
            )
            artists.append(r)
        return artists


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

        return (
            ax.pcolormesh(
                xedges,
                yedges,
                h.T,
                cmap=cmap,
                shading='auto'
            ),
            cmap
        )

    dh, dc = plot(decoy_points, cm.Reds)
    th, tc = plot(true_points, cm.Blues)

    # set norms (linear example)
    dh.set_norm(mpl.colors.Normalize(
        vmin=np.nanmin(dh.get_array()),
        vmax=np.nanmax(dh.get_array()),
    ))
    th.set_norm(mpl.colors.Normalize(
        vmin=np.nanmin(th.get_array()),
        vmax=np.nanmax(th.get_array()),
    ))

    decoy_sm = plt.cm.ScalarMappable(norm=dh.norm, cmap=dc)
    true_sm = plt.cm.ScalarMappable(norm=th.norm, cmap=tc)

    handles = [decoy_sm, true_sm]
    labels = ["Decoys", "True positives"]

    if long_points:
        long = ax.scatter(
            long_points.x, long_points.y,
            color=TOL_MAGENTA,
            marker='^',
        )
        handles.append(long)
        labels.append("Long sequences")

    ax.legend(
        handles=handles,
        labels=labels,
        handler_map={
            decoy_sm: HandlerColormap(dc, dh.norm),
            true_sm: HandlerColormap(tc, th.norm),
        },
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
