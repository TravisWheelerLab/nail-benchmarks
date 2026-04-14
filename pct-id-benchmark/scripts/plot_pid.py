#!/usr/bin/env python3
import argparse
from pathlib import Path

from plot import axes, Curve, COLORS, TOOL_COLORS, prefix_label, annotate, theta

import matplotlib.pyplot as plt
import numpy as np


X_MIN = 25
X_MAX = 10

Y_MIN = 0.0
Y_MAX = 1.0

X_LABEL = 17.5


MMSEQS_SEQ = [
    # "mmseqs-s5.7-ms2000.seq",
    # "mmseqs-s7.5-ms2000.seq",
    # "mmseqs-s10.0-ms2000.seq",
    "mmseqs-s12.0-ms2000.seq",
]

MMSEQS_PRF = [
    # "mmseqs-s5.7-ms2000.prf",
    # "mmseqs-s7.5-ms2000.prf",
    # "mmseqs-s10.0-ms2000.prf",
    "mmseqs-s12.0-ms2000.prf",
]

NAIL_SEQ = [
    # "nail-s5.7-ms2000.seq",
    # "nail-s7.5-ms2000.seq",
    # "nail-s10.0-ms2000.seq",
    "nail-s12.0-ms2000.seq",
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
    *NAIL_SEQ,
    *NAIL_PRF,
    *MMSEQS_SEQ,
    *MMSEQS_PRF,
]


def main(args):
    curves = []

    with open(args.pid) as f:
        lines = list(
            filter(lambda line: not line.startswith("#"), f.readlines()))

        fpr = float(lines[0].split()[-1])
        bin_cnt = Curve(lines[1])

        for line in lines[2:]:
            curve = Curve(line)
            if curve.prefix in PLOTTED:
                curve.auc = np.trapezoid(curve.y, curve.x)
                curves.append(curve)

    curves.sort(key=lambda c: c.auc, reverse=True)

    fig, ax = axes()

    plt.title("Recall by decreasing pairwise %identity")

    ax.set_zorder(2)
    ax.patch.set_visible(False)

    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)

    ax.set_xlabel("Decreasing pairwise %identity")
    ax.set_ylabel(f"Recall at {fpr} FP per search")

    ax.axhline(
        y=0.5,
        color="black",
        linestyle="--",
        linewidth=1, alpha=0.5
    )

    ax.grid(True)

    for curve in curves:
        (tool, _) = prefix_label(curve.prefix)

        color = TOOL_COLORS[curve.tool]

        ax.plot(
            curve.x, curve.y,
            color=color,
            marker='o',
            mfc='white',
            markersize=5,
        )

        def pos(x: float, curve: Curve):
            x1 = x + 1.0
            x2 = x - 1.0

            y = curve.approx_y(x)
            y1 = curve.approx_y(x1)
            y2 = curve.approx_y(x2)

            return ((x, y), theta(ax, x1, y1, x2, y2))

        pt = (X_LABEL, curve.approx_y(X_LABEL))
        offset = (0, 5)
        x = X_LABEL
        if curve.search_type == "prf":
            x = 16
        elif curve.search_type == "seq":
            if curve.tool == "hmmer":
                x = 16
            elif curve.tool == "nail":
                x = 17
                offset = (0, -20)
            elif curve.tool == "mmseqs":
                x = 14
                offset = (0, -20)
            elif curve.tool == "blast":
                x = 20
                offset = (0, -20)

        (pt, rotation) = pos(x, curve)
        annotate(
            ax, tool, pt, offset, color, rotation,
            linestyle='--',
            arrowstyle='-|>'
        )

    if args.bins:
        ax_bins = ax.twinx()
        ax_bins.set_zorder(1)
        ax_bins.set_ylabel("Sequence pairs (count)")

        ax_bins.bar(
            bin_cnt.x, bin_cnt.y,
            color=COLORS[-1],
            alpha=0.75,
            label="Sequence pair count in %ID bin"
        )

        fig.legend()

    plt.savefig(args.out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("pid", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--bins", action="store_true")
    args = p.parse_args()
    main(args)
