#!/usr/bin/env python3

from plot_common import parse_curve, plot_style

import bisect
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


X_MIN = 1e-3
X_MAX = 1_000.0
# X_MAX = 10.0

AUC_X_MAX = 1.0


def main(filename):
    curves = []

    with open(filename) as f:
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))

        for line in lines:
            label, x, y = parse_curve(line)
            z = bisect.bisect_right(x, X_MAX)
            x = x[:z]
            y = y[:z]

            zz = bisect.bisect_right(x, AUC_X_MAX)
            xx = x[:zz]
            yy = y[:zz]
            auc = np.trapezoid(yy, xx)

            curves.append((auc, label, x, y))

    curves.sort(key=lambda c: c[0], reverse=True)

    fig, ax = plt.subplots(figsize=(16, 9))
    plt.title("Recall by FPR")

    ax.set_xlabel("False positives per search (log scale)")
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_xscale("log")

    ax.set_ylabel("Recall")
    ax.set_ylim(0.0, 0.8)

    for _auc, label, x, y, in curves:
        ax.plot(
            x, y,
            **plot_style(label),
            markersize=3,
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
