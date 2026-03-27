#!/usr/bin/env python3
import argparse
from pathlib import Path

from plot import axes, Curve, COLORS, TOOL_COLORS, prefix_label

import matplotlib.pyplot as plt
import numpy as np


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
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))

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
    ax.set_xlim(25, 10)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("Decreasing pairwise %identity")
    ax.set_ylabel(f"Recall at {fpr} FP per search")

    ax.axhline(0.5, color="black", linestyle="--", linewidth=1)

    for curve in curves:
        tool, _params = prefix_label(curve.prefix)

        color = TOOL_COLORS[curve.tool]

        if "seq" in curve.extra:
            mfc = color
        else:
            mfc = 'white'

        ax.plot(
            curve.x, curve.y,
            color=color,
            marker='o',
            mfc=mfc,
            markersize=5,
            label=tool
        )

    # ax_bins = ax.twinx()
    # ax_bins.set_zorder(1)
    # ax_bins.set_ylabel("Sequence pairs (count)")

    # ax_bins.bar(
    #     bin_cnt.x, bin_cnt.y,
    #     color=COLORS[-1],
    #     alpha=0.75,
    #     label="pair count")

    fig.legend()

    plt.savefig(args.out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("pid", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    main(args)
