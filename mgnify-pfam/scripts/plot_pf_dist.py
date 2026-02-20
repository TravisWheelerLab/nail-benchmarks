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
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))

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

    return [s[-1], s[-2], s[0]]


def main(args):

    curves = parse(args.a)
    curves_b = parse(args.b)
    curves_c = parse(args.c)

    curves = list(filter(lambda c: c.pf_cnt > 1_000_000, curves))

    sample_a = sample_by_seed_cnt(curves)
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
        # ax.get_yaxis().set_major_formatter(lambda y, _: f"{y:g}")

    bins = np.arange(15, 100, 5)
    # bins = np.arange(15, 95, 10)

    labels = [
        "prefilter hits | -s 12.0",
        "prefilter hits | -s 10.0",
        "prefilter hits | -s 7.5",
    ]
    for (s_idx, curves) in enumerate(samples):
        for (q_idx, curve) in enumerate(curves):
            axes[q_idx].set_title(f"{queries[q_idx]}")
            axes[q_idx].set_xlabel("Prefilter score")
            axes[q_idx].set_ylabel("Number of hits")

            yy, edges = np.histogram(curve.data[psc], bins=bins, weights=curve.data[pf_cnt])

            ax = axes[q_idx]
            ax.plot(edges[:-1], yy,
                    color=colors[s_idx],
                    marker='o',
                    ms=4,
                    label=labels[s_idx])

            ax.set_xticks(edges[:-1])
            ax_labels = [
                f"{int(edges[i])}-{int(edges[i+1]-1)}"
                for i in range(len(edges) - 1)
            ]
            ax.set_xticklabels(ax_labels, rotation=45, fontsize=5)

    labels = [
        "nail seeds | -s 12.0",
        "nail seeds | -s 10.0",
        "nail seeds | -s 7.5",
    ]
    for (s_idx, curves) in enumerate(samples):
        for (q_idx, curve) in enumerate(curves):
            yy, xx = np.histogram(curve.data[psc], bins=bins, weights=curve.data[seed_cnt])
            axes[q_idx].plot(
                xx[:-1],
                yy,
                color=colors[s_idx],
                linestyle='--',
                marker='o',
                ms=4,
                markerfacecolor='white',
                label=labels[s_idx],
            )

    labels = [
        "nail hits | -s 12.0",
        "nail hits | -s 10.0",
        "nail hits | -s 7.5",
    ]
    for (s_idx, curves) in enumerate(samples):
        for (q_idx, curve) in enumerate(curves):
            yy, xx = np.histogram(curve.data[psc], bins=bins, weights=curve.data[nail_cnt])
            axes[q_idx].plot(
                xx[:-1],
                yy,
                color=colors[s_idx],
                linestyle=':',
                marker='d',
                ms=4,
                markerfacecolor='white',
                label=labels[s_idx],
            )

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
