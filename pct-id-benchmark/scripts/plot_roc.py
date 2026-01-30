#!/usr/bin/env python3

from plot import Curve, TOOL_COLORS

from pathlib import Path

import matplotlib.pyplot as plt


X_MIN = 1e-3
# X_MAX = 1_000.0
X_MAX = 1.0

AUC_X_MAX = 1.0


def main(filename):
    curves = []

    with open(filename) as f:
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))

        for line in lines:
            curves.append(Curve(line))

    fig, ax = plt.subplots(figsize=(16, 9))
    plt.title("Recall by FPR")

    ax.set_xlabel("False positives per search (log scale)")
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_xscale("log")

    ax.set_ylabel("Recall")
    ax.set_ylim(0.0, 0.8)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")

    ax.grid(True)

    for curve in curves:
        color = TOOL_COLORS[curve.tool]

        ax.plot(
            curve.x, curve.y,
            color=color,
        )

        for (xx, yy) in zip(curve.x, curve.y) :
            if xx >= X_MIN:
                y = yy
                break

        ax.text(
            X_MIN, y,
            s=curve.prefix,
            color=color,
            fontsize=8,
            fontweight="bold",
            ha="right",
            va="center",
        )

    plt.savefig(filename.with_suffix(".pdf"))


if __name__ == "__main__":
    import sys
    main(Path(sys.argv[1]))
