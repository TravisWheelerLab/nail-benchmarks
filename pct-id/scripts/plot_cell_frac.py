#!/usr/bin/env python3
import argparse
from pathlib import Path

from plot import axes, Scatter, TOL_CYAN

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


def heatmap(true_points, decoy_points, ax, long_points=None):
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

    dh = plot(decoy_points, cm.Reds)
    th = plot(true_points, cm.Blues)

    dh.set_norm(mpl.colors.Normalize(
        vmin=np.nanmin(dh.get_array()),
        vmax=np.nanmax(dh.get_array()),
    ))
    th.set_norm(mpl.colors.Normalize(
        vmin=np.nanmin(th.get_array()),
        vmax=np.nanmax(th.get_array()),
    ))

    decoy_sm = plt.cm.ScalarMappable(norm=dh.norm, cmap=dh.cmap)
    true_sm = plt.cm.ScalarMappable(norm=th.norm, cmap=th.cmap)

    handles = [decoy_sm, true_sm]
    labels = ["Decoys", "True positives"]

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
        handler_map={
            decoy_sm: HandlerColormap(dh.cmap, dh.norm),
            true_sm: HandlerColormap(th.cmap, th.norm),
        },
        handleheight=1.0,
        handlelength=6.0,
        labelspacing=1.2,
    )


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

    fig, ax = axes()

    heatmap(true_points, decoy_points, ax, long_points)

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
    p.add_argument("true_hits", type=Path)
    p.add_argument("decoy_hits", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--long_hits", type=Path)
    args = p.parse_args()
    main(args)
