#!/usr/bin/env python3
import argparse
from pathlib import Path

from plot import axes, Point, TOOL_COLORS, PREFIXES


import matplotlib.pyplot as plt

MMSEQS_SEQ = [
    "mmseqs-s5.7.seq",
    "mmseqs-s7.5.seq",
    "mmseqs-s10.0.seq",
    "mmseqs-s12.0.seq",
]

MMSEQS_PRF = [
    "mmseqs-s5.7.prf",
    "mmseqs-s7.5.prf",
    "mmseqs-s10.0.prf",
    "mmseqs-s12.0.prf",
]

NAIL = [
    "nail-ms2000.prf",
]

OTHER = [
    "blast.prf",
    "blast.seq",
    "hmmer.prf",
    "hmmer.seq",
]

PLOTTED = [
    *OTHER,
    *NAIL,
    *MMSEQS_SEQ,
    *MMSEQS_PRF,
]


def plot(args):
    points = []

    with open(args.time) as f:
        lines = list(filter(lambda line: not line.startswith("#"), f.readlines()))
        fpr = float(lines[0].split()[-1])

        for line in lines[1:]:
            p = Point(line)
            if p.prefix in PLOTTED:
                points.append(p)

    points.sort(key=lambda c: c.x, reverse=True)

    mmseqs_seq = list(filter(lambda p: p.prefix in MMSEQS_SEQ, points))
    mmseqs_prf = list(filter(lambda p: p.prefix in MMSEQS_PRF, points))

    ###

    fig, ax = axes()
    plt.title("Recall by runtime")

    ax.set_xlabel("Runtime (seconds; log scale)")
    ax.set_xscale("log")
    ax.set_xlim(1.0, 1e4)

    ax.set_ylabel(f"Recall at {fpr} FP per search")
    ax.set_ylim(0.2, 0.75)

    ax.grid(True)

    ###

    for p in points:
        color = TOOL_COLORS[p.tool]

        ax.scatter(
            p.x, p.y,
            color=color,
            linewidths=2,
            s=30,
        )

        (label, args) = PREFIXES[p.prefix]

        if p.prefix in OTHER:
            plt.text(
                p.x, p.y - 0.01,
                label,
                color=color,
                fontweight="bold",
                ha="left",
                va="top",
            )
    ###

    ##########################################

    def mmseqs_plot(pts, label):
        x = [p.x for p in pts]
        y = [p.y for p in pts]

        ax.plot(
            x, y,
            color=TOOL_COLORS['mmseqs'],
            linestyle='--'
        )

        plt.text(
            pts[-1].x, pts[-1].y - 0.01,
            label,
            color=TOOL_COLORS['mmseqs'],
            fontweight="bold",
            ha="center",
            va="top",
        )

        for p in pts:
            (_, args) = PREFIXES[p.prefix]
            plt.text(
                p.x, p.y + 0.01,
                f"{args}",
                color=TOOL_COLORS['mmseqs'],
                fontsize=10,
                ha="right",
                va="bottom",
            )

    ##########################################

    mmseqs_plot(mmseqs_seq, "mmseqs (seq)")
    mmseqs_plot(mmseqs_prf, "mmseqs (prf)")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("time", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    plot(args)
    plt.savefig(args.out)
