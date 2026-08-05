#!/usr/bin/env python3
import argparse
from pathlib import Path

from plot import axes, Scatter

import matplotlib.pyplot as plt


def main(args):
    with open(args.points) as f:
        points = Scatter(f.readlines())

    fig, ax = axes()

    ax.set_title("")
    ax.grid(True)

    ax.set_xlabel("Number of decoy hits (nail)")

    ax.set_ylabel("(GA - learned)")

    plt.scatter(
        points.x, points.y,
        alpha=0.25,
        edgecolors="none",
        linewidths=0,
        s=16,
    )

    plt.savefig(args.out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("points", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    main(args)
