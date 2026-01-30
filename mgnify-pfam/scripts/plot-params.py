#!/usr/bin/env python3

from plot import Point, link, partition

from pathlib import Path

import matplotlib.pyplot as plt


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
    import sys
    filename = Path(sys.argv[1])

    points = []

    with open(filename) as f:
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))

        for line in lines:
            points.append(Point(line))

    nail_points = list(filter(lambda p: p.params["tool"] == "nail", points))
    mmseqs_points = list(filter(lambda p: p.params["tool"] == "mmseqs", points))

    nail_points.sort(key=lambda p: p.x)
    mmseqs_points.sort(key=lambda p: p.x)

    # nail_links = link(nail_points, "s")
    # mmseqs_links = link(mmseqs_points, "s")

    fig, ax = plt.subplots(figsize=(16, 9))
    plt.title("")

    ax.set_xlabel("Runtime (seconds; log scale)")
    ax.set_xscale("log")
    ax.set_xlim(1.0, 1e4)

    ax.set_ylabel("")
    ax.set_ylim(0.0, 1.0)

    ax.grid(True)

    plot(nail_points, "blue")
    plot(mmseqs_points, "red")

    plt.savefig(filename.with_suffix(".pdf"))
