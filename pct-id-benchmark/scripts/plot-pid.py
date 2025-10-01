#!/usr/bin/env python3

from plot_common import parse_curve, plot_style, COLORS

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main(filename):
    curves = []

    with open(filename) as f:
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))

        fpr = float(lines[0].split()[-1])
        bin_cnt = parse_curve(lines[1])

        for line in lines[2:]:
            label, x, y = parse_curve(line)
            auc = np.trapezoid(y, x)
            curves.append((auc, label, x, y))

    curves.sort(key=lambda c: c[0], reverse=True)

    fig, ax_recall = plt.subplots(figsize=(16, 9))
    plt.title("Recall by decreasing pairwise %identity")

    ax_recall.set_zorder(2)
    ax_recall.patch.set_visible(False)
    ax_recall.set_xlim(25, 10)
    ax_recall.set_ylim(0.0, 1.0)
    ax_recall.set_xlabel("Decreasing pairwise %identity")
    ax_recall.set_ylabel(f"Recall at {fpr} FP per search")

    ax_recall.axhline(0.5, color="black", linestyle="--", linewidth=1)

    for _auc, label, x, y, in curves:
        ax_recall.plot(
            x, y,
            **plot_style(label),
            markersize=5,
        )

    ax_bins = ax_recall.twinx()
    ax_bins.set_zorder(1)
    ax_bins.set_ylabel("sequence pairs (count)")

    ax_bins.bar(bin_cnt[1], bin_cnt[2], color=COLORS[-1], alpha=0.75, label="pair count")

    fig.legend()

    plt.savefig(filename.with_suffix(".pdf"))


if __name__ == "__main__":
    import sys
    main(Path(sys.argv[1]))
