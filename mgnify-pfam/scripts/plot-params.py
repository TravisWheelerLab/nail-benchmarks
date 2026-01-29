#!/usr/bin/env python3

from pathlib import Path
import re

import matplotlib.pyplot as plt


class Point:
    x: float
    y: float
    label: str
    params: dict[str, float]

    def __init__(self, line: str):
        prefix, search_type, point = line.split(",", 2)

        x, y = map(float, re.match(r'\(\s*([-+]?\d*\.?\d+)\s*,\s*([-+]?\d*\.?\d+)\s*\)', point).groups())
        self.x = float(x)
        self.y = float(y)

        tokens = prefix.split('-')

        self.params = {}
        self.params["tool"] = tokens[0]

        for tok in tokens[1:]:
            m = re.match(r'([a-zA-Z]+)([0-9.]+)', tok)
            p, v = m.groups()
            self.params[p] = float(v)


def link(points: [Point], key: str):
    aa = list(set([p.params[key] for p in points]))
    aaa = [
        [
            [p.x for p in points if p.params[key] == v],
            [p.y for p in points if p.params[key] == v],
            f"{key}-{v}"
        ]
        for v in aa
    ]

    return aaa


def main(filename):
    points = []

    with open(filename) as f:
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))

        for line in lines:
            points.append(Point(line))

    nail_points = list(filter(lambda p: p.params["tool"] == "nail", points))
    mmseqs_points = list(filter(lambda p: p.params["tool"] == "mmseqs", points))

    nail_points.sort(key=lambda p: p.x)
    mmseqs_points.sort(key=lambda p: p.x)

    nail_links = link(nail_points, "s")
    mmseqs_links = link(mmseqs_points, "s")

    fig, ax = plt.subplots(figsize=(16, 9))
    plt.title("")

    ax.set_xlabel("Runtime (seconds; log scale)")
    ax.set_xscale("log")
    ax.set_xlim(1.0, 1e4)

    ax.set_ylabel("")
    ax.set_ylim(0.0, 1.0)

    ax.grid(True)

    for (x, y, l) in nail_links:
        ax.plot(x, y, color="blue")
        plt.text(
            x[0], y[0],
            l,
            color='blue',
            fontsize=8,
            fontweight="bold",
            ha="right",
            va="bottom",
        )

    for p in nail_points:
        ax.scatter(
            p.x, p.y,
            color="blue",
            linewidths=2,
            s=30,
        )

    for (x, y, l) in mmseqs_links:
        ax.plot(x, y, color="red")
        plt.text(
            x[0], y[0],
            l,
            color='red',
            fontsize=8,
            fontweight="bold",
            ha="right",
            va="bottom",
        )

    for p in mmseqs_points:
        ax.scatter(
            p.x, p.y,
            color="red",
            linewidths=2,
            s=30,
        )

    plt.savefig(filename.with_suffix(".pdf"))


if __name__ == "__main__":
    import sys
    main(Path(sys.argv[1]))
