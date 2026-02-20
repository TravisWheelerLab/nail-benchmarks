#!/usr/bin/env python3
import argparse

from plot import Scatter, TOL_CYAN, TOL_MAGENTA

from pathlib import Path

import matplotlib.pyplot as plt


def main(args):

    with open(args.true_hits) as f:
        true_points = Scatter(f.readlines())

    with open(args.decoy_hits) as f:
        decoy_points = Scatter(f.readlines())

    fig, ax = plt.subplots(figsize=(16, 9))
    plt.title("Cells computed")

    ax.set_xlabel("Total DP cells (product of target and query lengths)")
    ax.set_xscale("log")
    ax.set_xlim(1e2, 1e10)

    ax.set_ylabel("Fraction of DP cells computed by sparse F/B")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1.1)

    ax.grid(True)

    alpha = 0.6

    ax.scatter(
        decoy_points.x, decoy_points.y,
        label="decoys",
        color=TOL_MAGENTA,
        marker='v',
        alpha=alpha,
        rasterized=True,
    )

    ax.scatter(
        true_points.x, true_points.y,
        label="True positives",
        color=TOL_CYAN,
        marker='^',
        alpha=alpha,
        rasterized=True,
    )

    plt.legend()
    plt.savefig(args.true_hits.with_name("cells.pdf"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("true_hits", type=Path)
    p.add_argument("decoy_hits", type=Path)
    args = p.parse_args()
    main(args)
