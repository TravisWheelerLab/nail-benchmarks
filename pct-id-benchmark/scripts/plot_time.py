#!/usr/bin/env python3
import argparse
from pathlib import Path

from plot import axes, Point, TOOL_COLORS, prefix_label


import matplotlib.pyplot as plt
import numpy as np


MMSEQS_SEQ = [
    "mmseqs-s5.7-ms2000.seq",
    "mmseqs-s7.5-ms2000.seq",
    "mmseqs-s10.0-ms2000.seq",
    "mmseqs-s12.0-ms2000.seq",
    # "mmseqs-s14.0-ms2000.seq",
]

MMSEQS_PRF = [
    "mmseqs-s5.7-ms2000.prf",
    "mmseqs-s7.5-ms2000.prf",
    "mmseqs-s10.0-ms2000.prf",
    "mmseqs-s12.0-ms2000.prf",
    "mmseqs-s14.0-ms2000.prf",
]

NAIL_SEQ = [
    "nail-s5.7-ms2000.seq",
    "nail-s7.5-ms2000.seq",
    "nail-s10.0-ms2000.seq",
    "nail-s12.0-ms2000.seq",
    # "nail-s14.0-ms2000.seq",
]

NAIL_PRF = [
    "nail-s5.7-ms2000.prf",
    "nail-s7.5-ms2000.prf",
    "nail-s10.0-ms2000.prf",
    "nail-s12.0-ms2000.prf",
    "nail-s14.0-ms2000.prf",
]

OTHER = [
    "blast.prf",
    "blast.seq",
    "hmmer.prf",
    "hmmer.seq",
]

PLOTTED = [
    *OTHER,
    *NAIL_SEQ,
    *NAIL_PRF,
    *MMSEQS_SEQ,
    *MMSEQS_PRF,
]

ax = None


def plot_line(pts, label, layout=None, label_layout=None):
    global ax

    color = TOOL_COLORS[label.split()[0]]

    x = np.array([p.x for p in pts])
    y = np.array([p.y for p in pts])

    ax.plot(
        x, y,
        color=color,
        linestyle='--'
    )

    if label_layout is None:
        o = (10, -10)
        ha = "center"
    else:
        o = label_layout[0]
        ha = label_layout[1]

    last_pt = (pts[-1].x, pts[-1].y)
    ax.annotate(
        label,
        last_pt,
        xytext=o,
        textcoords="offset pixels",
        color=color,
        fontweight="bold",
        ha=ha,
    )

    if layout is None:
        offsets = [(10, -10) for _ in range(len(pts))]
        ha_list = ["left" for _ in range(len(pts))]
    else:
        offsets = [l[0] for l in layout]
        ha_list = [l[1] for l in layout]

    for (p, o, ha) in zip(pts, offsets, ha_list):
        _, args = prefix_label(p.prefix, exclude=["ms"], reformat=[
                               ("--mmseqs-s", "-s")])
        args = " ".join(args)

        pt = (p.x, p.y)
        ax.annotate(
            f"{args}",
            pt,
            xytext=o,
            textcoords="offset pixels",
            color=color,
            fontsize=10,
            ha=ha,
        )


def plot(args):
    points = []

    with open(args.time) as f:
        lines = list(
            filter(lambda line: not line.startswith("#"), f.readlines()))
        fpr = float(lines[0].split()[-1])

        for line in lines[1:]:
            p = Point(line)
            if p.prefix in PLOTTED:
                points.append(p)

    points.sort(key=lambda c: c.x, reverse=True)

    nail_prf = list(filter(lambda p: p.prefix in NAIL_PRF, points))
    nail_seq = list(filter(lambda p: p.prefix in NAIL_SEQ, points))
    mmseqs_seq = list(filter(lambda p: p.prefix in MMSEQS_SEQ, points))
    mmseqs_prf = list(filter(lambda p: p.prefix in MMSEQS_PRF, points))

    ###

    global ax, fig
    fig, ax = axes()
    plt.title("Recall by runtime")

    ax.set_xlabel("Runtime (seconds; log scale)")
    ax.set_xscale("log")

    ax.set_ylabel(f"Recall at {fpr} FP per search")

    ax.set_xlim(1.0, 1e3)
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

        (label, args) = prefix_label(p.prefix, exclude=["ms"])
        pt = (p.x, p.y)

        if p.prefix == "blast.prf":
            ax.annotate(
                label,
                pt,
                xytext=(-30, 15),
                textcoords="offset pixels",
                color=color,
                fontweight="bold",
                ha="center",
            )
        elif p.prefix == "blast.seq":
            ax.annotate(
                label,
                pt,
                xytext=(-10, -50),
                textcoords="offset pixels",
                color=color,
                fontweight="bold",
                ha="center",
                arrowprops=dict(
                    arrowstyle="->",
                    color=color,
                )
            )
        elif p.prefix == "hmmer.prf":
            ax.annotate(
                label,
                pt,
                xytext=(0, 15),
                textcoords="offset pixels",
                color=color,
                fontweight="bold",
                ha="center",
            )
        elif p.prefix == "hmmer.seq":
            ax.annotate(
                label,
                pt,
                xytext=(-75, 15),
                textcoords="offset pixels",
                color=color,
                fontweight="bold",
                ha="center",
            )

    plot_line(
        nail_prf,
        "nail (profile)",
        layout=[
            [(0, 12), "center"],
            [(0, 12), "center"],
            [(10, -10), "left"],
            [(10, -2.5), "left"],
            [(10, -2.5), "left"],
        ],
        label_layout=[(-50, 260), "left"]
    )

    plot_line(
        nail_seq,
        "nail (sequence)",
        layout=[
            [(0, -15), "right"],
            [(0, -25), "center"],
            [(10, -10), "left"],
            [(10, -2.5), "left"],
        ],
        label_layout=[(175, -30), "right"]
    )

    plot_line(
        mmseqs_prf,
        "mmseqs (profile)",
        layout=[
            [(0, -17.5), "center"],
            [(0, -17.5), "center"],
            [(0, 10), "center"],
            [(10, -10), "left"],
            [(10, 0), "left"],
        ],
        label_layout=[(-10, 50), "right"]
    )

    plot_line(
        mmseqs_seq,
        "mmseqs (sequence)",
        layout=[
            [(0, 10), "center"],
            [(0, 10), "center"],
            [(-10, 0), "right"],
            [(-10, 0), "right"],
        ],
        label_layout=[(-210, -30), "left"]
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("time", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    plot(args)
    plt.savefig(args.out)
