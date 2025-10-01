#!/usr/bin/env python3

from plot_common import parse_curve, plot_style

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main(filename):
    curves = []

    with open(filename) as f:
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))

        for line in lines:
            label, x, y = parse_curve(line)
            auc = np.trapezoid(y, x)
            curves.append((auc, label, x, y))

    curves.sort(key=lambda c: c[0], reverse=True)

    fig, ax = plt.subplots(figsize=(16, 9))
    plt.title("Recall by decreasing pairwise %identity")

    ax.set_xlim(1e-3, 5.0)
    ax.set_xscale("log")
    ax.set_xlabel("")
    ax.set_ylabel("")

    for _auc, label, x, y, in curves:
        ax.plot(
            x, y,
            **plot_style(label),
            markersize=3,
        )

    fig.legend()

    plt.savefig(filename.with_suffix(".pdf"))


if __name__ == "__main__":
    import sys
    main(Path(sys.argv[1]))
