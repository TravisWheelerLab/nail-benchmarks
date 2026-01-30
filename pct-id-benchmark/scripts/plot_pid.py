#!/usr/bin/env python3

from plot import Curve, COLORS, TOOL_COLORS

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main(filename):
    curves = []

    with open(filename) as f:
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))

        fpr = float(lines[0].split()[-1])
        bin_cnt = Curve(lines[1])

        for line in lines[2:]:
            curve = Curve(line)
            curve.auc = np.trapezoid(curve.y, curve.x)
            curves.append(curve)

    curves.sort(key=lambda c: c.auc, reverse=True)

    fig, ax_recall = plt.subplots(figsize=(16, 9))
    plt.title("Recall by decreasing pairwise %identity")

    ax_recall.set_zorder(2)
    ax_recall.patch.set_visible(False)
    ax_recall.set_xlim(25, 10)
    ax_recall.set_ylim(0.0, 1.0)
    ax_recall.set_xlabel("Decreasing pairwise %identity")
    ax_recall.set_ylabel(f"Recall at {fpr} FP per search")

    ax_recall.axhline(0.5, color="black", linestyle="--", linewidth=1)

    for curve in curves:
        color = TOOL_COLORS[curve.tool]

        if "seq" in curve.extra:
            mfc = color
            label = f"{curve.prefix} seq"
        elif "prf" in curve.extra:
            mfc = 'white'
            label = f"{curve.prefix} prf"

        ax_recall.plot(
            curve.x, curve.y,
            color=color,
            marker='o',
            mfc=mfc,
            markersize=5,
            label=label
        )

    ax_bins = ax_recall.twinx()
    ax_bins.set_zorder(1)
    ax_bins.set_ylabel("sequence pairs (count)")

    ax_bins.bar(bin_cnt.x, bin_cnt.y, color=COLORS[-1], alpha=0.75, label="pair count")

    fig.legend(
        fontsize=8,
        markerscale=1.0,
    )

    plt.savefig(filename.with_suffix(".pdf"))


if __name__ == "__main__":
    import sys
    main(Path(sys.argv[1]))
