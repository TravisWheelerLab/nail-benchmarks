#!/usr/bin/env python3

from plot import Point, link, partition, TOOL_COLORS
from plot_time import plot as plot_other

import argparse
from pathlib import Path

import matplotlib.pyplot as plt


def parse(filename):
    points = []
    with open(filename) as f:
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))
        fpr = float(lines[0].split()[-1])

        for line in lines[1:]:
            points.append(Point(line))

    return (points, fpr)


def plot(points, color):
    (prog_points, ms_points) = partition(points, lambda p: "prog" in p.extra)
    links = link(points, "s")

    for (x, y, l) in links:
        ax.plot(
            x, y,
            color=color,
            linestyle='--',
            linewidth=0.5,
        )

        plt.text(
            x[0] + 0.01, y[0] + 0.01,
            l,
            color=color,
            fontsize=6,
            fontweight="bold",
            ha="right",
            va="top",
        )

    for p in ms_points:
        ax.scatter(
            p.x, p.y,
            color=color,
            s=3,
        )

    for p in prog_points:
        ax.scatter(
            p.x, p.y,
            color=color,
            marker="^",
            s=12,
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("nail", type=Path)
    p.add_argument("mmseqs", type=Path)
    p.add_argument('--overlay', default=None)
    p.add_argument("out", type=Path)

    args = p.parse_args()

    nail_points, nail_fpr = parse(args.nail)
    mmseqs_points, mmseqs_fpr = parse(args.mmseqs)

    nail_points.sort(key=lambda p: p.x)
    mmseqs_points.sort(key=lambda p: p.x)

    assert (nail_fpr == mmseqs_fpr)
    fpr = nail_fpr

    fig, ax = plt.subplots(figsize=(16, 9))
    plt.title("Recall by runtime")

    ax.set_xlabel("Runtime (seconds; log scale)")
    ax.set_xscale("log")
    ax.set_xlim(1.0, 1e4)

    ax.set_ylabel(f"Recall at {fpr} FP per search")
    ax.set_ylim(0.0, 0.8)
    xp = 1.5
    ax.set_yscale('function', functions=(lambda x: x**xp, lambda x: x**(1 / xp)))
    ax.grid(True)

    if args.overlay is not None:
        plot_other(
            args.overlay, ax,
            exclude=[lambda p: "seq" in p.extra]
        )

    plot(nail_points, TOOL_COLORS["nail"])
    plot(mmseqs_points, TOOL_COLORS["mmseqs"])

    plt.savefig(args.out)
