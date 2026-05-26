#!/usr/bin/env python3

from plot import CurveND, COLORS

from pathlib import Path
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np

psc = 0
seed_frac = 1
seed_cnt = 2
nail_frac = 3
nail_cnt = 4
pf_cnt = 5


def parse(filename):
    curves = []
    with open(filename) as f:
        lines = list(
            filter(lambda line: not line.startswith("#"), f.readlines()))

        for line in lines:
            curve = CurveND(line, 6)
            curve.seed_cnt = sum(curve.data[seed_cnt])
            curve.pf_cnt = sum(curve.data[pf_cnt])
            curves.append(curve)

    return curves


def sample_by_seed_cnt(curves, k=15, w=1):
    curves.sort(key=lambda c: c.seed_cnt)
    sz = len(curves)
    idx = (np.linspace(0, 1, k) ** w * (sz - 1)).astype(int)
    s = [curves[i] for i in idx]

    return [s[-2], s[-1], s[0]]


def main(args):

    curves_a = parse(args.a)
    curves_b = parse(args.b)
    curves_c = parse(args.c)

    curves_a = list(filter(lambda c: c.pf_cnt > 1_000_000, curves_a))

    sample_a = sample_by_seed_cnt(curves_a)
    queries = [c.prefix for c in sample_a]

    sample_b = [c for c in curves_b if c.prefix in queries]
    sample_c = [c for c in curves_c if c.prefix in queries]

    sample_a.sort(key=lambda c: queries.index(c.prefix))
    sample_b.sort(key=lambda c: queries.index(c.prefix))
    sample_c.sort(key=lambda c: queries.index(c.prefix))

    samples = [sample_a, sample_b, sample_c]

    ##########################################################################

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    colors = [COLORS[0], COLORS[1], COLORS[5]]

    for ax in axes.flat:
        ax.set_xlim(15, 90)
        ax.set_ylim(0, 1e6)
        ax.set_yscale("symlog", linthresh=1)
        ax.set_yticks([0, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6])

    bins = np.arange(15, 100, 5)
    sens_labels = ["-s 12.0", "-s 10.0", "-s 7.5"]

    series_cfg = [
        {
            "col": pf_cnt,
            "prefix": "prefilter hits",
            "linestyle": "-",
            "marker": "o",
        },
        {
            "col": seed_cnt,
            "prefix": "nail seeds",
            "linestyle": "--",
            "marker": "o",
            "markerfacecolor": "white",
        },
        {
            "col": nail_cnt,
            "prefix": "nail hits",
            "linestyle": ":",
            "marker": "d",
            "markerfacecolor": "white",
        },
    ]

    for cfg in series_cfg:
        plot_kw = {k: v for k, v in cfg.items() if k not in ("col", "prefix")}
        for s_idx, curves in enumerate(samples):
            for q_idx, curve in enumerate(curves):
                yy, _ = np.histogram(
                    curve.data[psc],
                    bins=bins,
                    weights=curve.data[cfg["col"]]
                )

                label = f"{cfg['prefix']} | {sens_labels[s_idx]}"

                axes[q_idx].plot(
                    bins[:-1],
                    yy,
                    color=colors[s_idx],
                    ms=4,
                    label=label,
                    **plot_kw,
                )

    ax_labels = [f"{int(bins[i])}-{int(bins[i+1]-1)
                                   }" for i in range(len(bins) - 1)]
    for q_idx, query in enumerate(queries):
        ax = axes[q_idx]
        ax.set_title(query)
        ax.set_xlabel("Prefilter score")
        ax.set_ylabel("Number of hits")
        ax.set_xticks(bins[:-1])
        ax.set_xticklabels(ax_labels, rotation=45, fontsize=5)

    plt.legend()

    ##########################################################################

    plt.savefig("prefilter_scores.pdf")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("a", type=Path)
    p.add_argument("b", type=Path)
    p.add_argument("c", type=Path)

    args = p.parse_args()
    main(args)
