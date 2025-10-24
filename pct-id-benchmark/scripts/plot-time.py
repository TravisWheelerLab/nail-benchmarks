#!/usr/bin/env python3

from plot_common import parse_point, scatter_style

from pathlib import Path

import matplotlib.pyplot as plt


def main(filename):
    points = []

    with open(filename) as f:
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))
        fpr = float(lines[0].split()[-1])

        for line in lines[1:]:
            label, x, y = parse_point(line)
            points.append((label, x, y))

    points.sort(key=lambda c: c[1], reverse=True)

    fig, ax = plt.subplots(figsize=(16, 9))
    plt.title("Recall by runtime")

    ax.set_xlabel("Runtime (seconds; log scale)")
    ax.set_xscale("log")
    ax.set_xlim(1.0, 1e4)

    ax.set_ylabel(f"Recall at {fpr} FP per search")
    ax.set_ylim(0.0, 0.8)

    ax.grid(True)

    for label, x, y, in points:
        ax.scatter(
            x, y,
            **scatter_style(label),
            s=20,
        )

    fig.legend(
        fontsize=8,
        markerscale=1.0,
        loc='lower right'
    )

    plt.savefig(filename.with_suffix(".pdf"))


if __name__ == "__main__":
    import sys
    main(Path(sys.argv[1]))
