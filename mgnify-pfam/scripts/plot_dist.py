#!/usr/bin/env python3
import argparse
from pathlib import Path

from plot import axes, Point

import matplotlib.pyplot as plt
import numpy as np

def main(args):
    points = []
    with open(args.points) as f:
        for line in f:
            point = Point(line)
            points.append(point)

    fig, ax = axes()

    ax.set_title("")
    ax.grid(True)

    ax.set_xlabel("Number of decoy hits")
    ax.set_xscale("log")

    ax.set_ylabel("Query count")

    nail = [p.x for p in points]
    mmseqs = [p.y for p in points]

    min_val = min(*nail, *mmseqs)
    max_val = max(*nail, *mmseqs)

    bins = np.arange(min_val, max_val + 1, 1)

    plt.hist(
        mmseqs,
        bins=bins,
        alpha=0.50,
        label="mmseqs"
    )

    plt.hist(
        nail,
        bins=bins,
        alpha=0.50,
        label="nail"
    )

    ax.legend()

    plt.savefig(args.out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("points", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    main(args)
