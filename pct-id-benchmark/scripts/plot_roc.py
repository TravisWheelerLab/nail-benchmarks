#!/usr/bin/env python3
import argparse
from pathlib import Path

from plot import axes, Curve, TOOL_COLORS, PREFIXES


import matplotlib.pyplot as plt


X_MIN = 1e-3
# X_MAX = 1_000.0
X_MAX = 1.0

AUC_X_MAX = 1.0

MMSEQS_SEQ = [
    # "mmseqs-s5.7.seq",
    # "mmseqs-s7.5.seq",
    # "mmseqs-s10.0.seq",
    "mmseqs-s12.0.seq",
]

MMSEQS_PRF = [
    # "mmseqs-s5.7.prf",
    # "mmseqs-s7.5.prf",
    # "mmseqs-s10.0.prf",
    "mmseqs-s12.0.prf",
]

NAIL_PRF = [
    # "nail-s5.7-ms2000.prf",
    # "nail-s7.5-ms2000.prf",
    # "nail-s10.0-ms2000.prf",
    "nail-s12.0-ms2000.prf",
]

OTHER = [
    "blast.prf",
    "blast.seq",
    "hmmer.prf",
    "hmmer.seq",
]

PLOTTED = [
    *OTHER,
    *NAIL_PRF,
    *MMSEQS_SEQ,
    *MMSEQS_PRF,
]


def main(args):
    curves = []

    with open(args.roc) as f:
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))

        for line in lines:
            curve = Curve(line)
            if curve.prefix in PLOTTED:
                curves.append(curve)

    fig, ax = axes()

    plt.title("Recall by FPR")

    ax.set_xlabel("False positives per search (log scale)")
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_xscale("log")

    ax.set_ylabel("Recall")
    ax.set_ylim(0.35, 0.8)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")

    ax.grid(True)

    for curve in curves:
        color = TOOL_COLORS[curve.tool]
        (label, params) = PREFIXES[curve.prefix]

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
            s=label,
            color=color,
            fontsize=8,
            fontweight="bold",
            ha="right",
            va="center",
        )

    plt.savefig(args.out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("roc", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    main(args)
