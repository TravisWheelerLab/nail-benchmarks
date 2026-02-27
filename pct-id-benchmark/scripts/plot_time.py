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

NAIL_PRF = [
    # "nail-s5.7-ms2000.prf",
    "nail-s7.5-ms2000.prf",
    "nail-s10.0-ms2000.prf",
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
    *NAIL_PRF,
    *MMSEQS_SEQ,
    *MMSEQS_PRF,
]

ax = None


def nail_plot(pts, label):
    global ax

    color = TOOL_COLORS['nail']
    x = [p.x for p in pts]
    y = [p.y for p in pts]

    ax.plot(
        x, y,
        color=color,
        linestyle='--'
    )

    plt.text(
        pts[-1].x, pts[-1].y - 0.01,
        label,
        color=color,
        fontweight="bold",
        ha="center",
        va="top",
    )

    for p in pts:
        (_, args) = PREFIXES[p.prefix]
        plt.text(
            p.x, p.y + 0.005,
            f"{args}",
            color=color,
            fontsize=10,
            ha="right",
            va="bottom",
        )


def mmseqs_plot(pts, label):
    global ax
    x = [p.x for p in pts]
    y = [p.y for p in pts]

    color = TOOL_COLORS['mmseqs']

    ax.plot(
        x, y,
        color=color,
        linestyle='--'
    )

    plt.text(
        pts[-1].x, pts[-1].y - 0.01,
        label,
        color=color,
        fontweight="bold",
        ha="center",
        va="top",
    )

    for p in pts:
        (_, args) = PREFIXES[p.prefix]
        plt.text(
            p.x, p.y + 0.005,
            f"{args}",
            color=color,
            fontsize=10,
            ha="right",
            va="bottom",
        )


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

    nail_prf = list(filter(lambda p: p.prefix in NAIL_PRF, points))
    mmseqs_seq = list(filter(lambda p: p.prefix in MMSEQS_SEQ, points))
    mmseqs_prf = list(filter(lambda p: p.prefix in MMSEQS_PRF, points))

    ###

    global ax
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
                ha="center",
                va="top",
            )
    ###

    nail_plot(nail_prf, "nail (prf)")
    mmseqs_plot(mmseqs_seq, "mmseqs (seq)")
    mmseqs_plot(mmseqs_prf, "mmseqs (prf)")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("time", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    plot(args)
    plt.savefig(args.out)
