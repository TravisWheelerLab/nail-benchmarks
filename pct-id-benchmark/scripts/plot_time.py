#!/usr/bin/env python3
import argparse
from pathlib import Path

from plot import axes, Point, TOOL_COLORS


import matplotlib.pyplot as plt


def plot(args, ax=None, exclude=[]):
    points = []

    with open(args.time) as f:
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))
        fpr = float(lines[0].split()[-1])

        for line in lines[1:]:
            points.append(Point(line))

    points.sort(key=lambda c: c.x, reverse=True)

    if ax is None:
        fig, ax = axes()
        plt.title("Recall by runtime")

        ax.set_xlabel("Runtime (seconds; log scale)")
        ax.set_xscale("log")
        ax.set_xlim(1.0, 1e4)

        ax.set_ylabel(f"Recall at {fpr} FP per search")
        ax.set_ylim(0.0, 0.8)

        ax.grid(True)


    for point in points:
        if any(e(point) for e in exclude):
            continue

        color = TOOL_COLORS[point.tool]

        ax.scatter(
            point.x, point.y,
            color=color,
            linewidths=2,
            s=30,
        )

        plt.text(
            point.x, point.y + 0.01,
            point.prefix,
            color=color,
            fontweight="bold",
            ha="center",
            va="bottom",
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("time", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    plot(args)
    plt.savefig(args.out)
