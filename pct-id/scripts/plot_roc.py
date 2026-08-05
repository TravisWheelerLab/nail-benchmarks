#!/usr/bin/env python3
import argparse
from pathlib import Path

from plot import axes, Curve, TOOL_COLORS, prefix_label, annotate


import matplotlib.pyplot as plt


X_MIN = 1e-3
X_MAX = 1.0

Y_MIN = 0.35
Y_MAX = 0.8

X_LABEL = 3.0e-2

AUC_X_MAX = 1.0

MMSEQS_SEQ = [
    # "mmseqs-s5.7-ms2000.seq",
    "mmseqs-s7.5-ms2000.seq",
    # "mmseqs-s10.0-ms2000.seq",
    "mmseqs-s12.0-ms2000.seq",
]

MMSEQS_PRF = [
    # "mmseqs-s5.7-ms2000.prf",
    "mmseqs-s7.5-ms2000.prf",
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
    "diamond-ultra-sens.seq",
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

    with open(args.roc) as f:
        lines = list(
            filter(lambda line: not line.startswith("#"), f.readlines()))

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
    ax.set_ylim(Y_MIN, Y_MAX)

    ax.grid(True)

    for curve in curves:
        color = TOOL_COLORS[curve.tool]
        tool, _params = prefix_label(curve.prefix)

        ax.plot(
            curve.x, curve.y,
            color=color,
            drawstyle="steps-post",
        )

        # a curve that never reaches the annotation x cannot be labelled; draw
        # it unlabelled rather than losing the whole figure to one short curve
        try:
            pt = (X_LABEL, curve.approx_y(X_LABEL))
            offset = (0, 0)
            rotation = 0
            if curve.search_type == "prf":
                if curve.tool == "hmmer":
                    offset = (0, 20)
                elif curve.tool == "nail":
                    offset = (0, 20)
                elif curve.tool == "mmseqs":
                    offset = (0, -25)
                elif curve.tool == "blast":
                    offset = (0, 20)
            elif curve.search_type == "seq":
                if curve.tool == "hmmer":
                    offset = (0, 20)
                elif curve.tool == "nail":
                    x = X_LABEL * 2.0
                    pt = (x, curve.approx_y(x))
                    offset = (0, -30)
                elif curve.tool == "mmseqs":
                    x = X_LABEL / 2.0
                    pt = (x, curve.approx_y(x))
                    offset = (0, -30)
                elif curve.tool == "blast":
                    offset = (0, 15)
                    rotation = 10
                elif curve.tool == "diamond":
                    x = X_LABEL * 10.0
                    pt = (x, curve.approx_y(x))
                    offset = (0, -30)

            annotate(
                ax, tool, pt, offset, color, rotation,
                linestyle='--',
                arrowstyle='-|>'
            )
        except ValueError:
            pass

    plt.savefig(args.out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("roc", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    main(args)
