#!/usr/bin/env python3

from plot_common import parse_point

import argparse
from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np


def parse(filename):
    points = []
    with open(filename) as f:
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))
        fpr = float(lines[0].split()[-1])

        for line in lines[1:]:
            label, x, y = parse_point(line)
            s = float(label.split('-')[1][1:])
            points.append((label, s, x, y))

    return (points, fpr)


def plot(points, c):
    s_list = list(set([s for (_, s, _, _) in points]))
    s_list.sort()
    colors = c(np.linspace(0.3, 1.0, len(s_list)))
    colors = {l: c for (l, c) in zip(s_list, colors)}

    groups = defaultdict(list)
    for p in points:
        groups[p[1]].append(p)

    for s in s_list:
        g = groups[s]
        color = colors[s]
        g.sort(key=lambda x: x[2])

        x = [x for *_, x, _ in g]
        y = [y for *_, _, y in g]

        plt.plot(
            x, y,
            color=color,
            marker='o',
            markersize='3'
        )

        plt.text(
            x[0], y[0],
            f"s{s}",
            color=color,
            fontsize=8,
            fontweight="bold",
            ha="right",
            va="bottom",
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("nail", type=Path)
    p.add_argument("mmseqs", type=Path)
    p.add_argument("out", type=Path)

    args = p.parse_args()

    nail_points, nail_fpr = parse(args.nail)
    mmseqs_points, mmseqs_fpr = parse(args.mmseqs)

    assert (nail_fpr == mmseqs_fpr)

    fpr = nail_fpr

    fig, ax = plt.subplots(figsize=(16, 9))
    plt.title("Recall by runtime")

    ax.set_xlabel("Runtime (seconds; log scale)")
    ax.set_xscale("log")
    ax.set_xlim(1.0, 1e4)

    ax.set_ylabel(f"Recall at {fpr} FP per search")
    ax.set_ylim(0.0, 0.8)

    ax.grid(True)

    plot(nail_points, plt.cm.Blues)
    plot(mmseqs_points, plt.cm.Reds)

    plt.savefig(args.out)
