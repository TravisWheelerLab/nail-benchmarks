#!/usr/bin/env python3

from plot import Point, axes, TOOL_COLORS, prefix_label, annotate

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# NAIL_S57 = [
#     "nail-s5.7-prog.prf",
#     "nail-s5.7-ms300.prf",
#     "nail-s5.7-ms2000.prf",
#     "nail-s5.7-ms2147483647.prf",
# ]

# NAIL_S75 = [
#     "nail-s7.5-prog.prf",
#     "nail-s7.5-ms300.prf",
#     "nail-s7.5-ms2000.prf",
#     "nail-s7.5-ms2147483647.prf",
# ]

# NAIL_S90 = [
#     "nail-s9.0-prog.prf",
#     "nail-s9.0-ms300.prf",
#     "nail-s9.0-ms300.prf",
#     "nail-s9.0-ms2000.prf",
# ]

# NAIL_S100 = [
#     "nail-s10.0-prog.prf",
#     "nail-s10.0-ms300.prf",
#     "nail-s10.0-ms2000.prf",
#     "nail-s10.0-ms2147483647.prf",
# ]

# NAIL_S110 = [
#     "nail-s11.0-prog.prf",
#     "nail-s11.0-ms300.prf",
#     "nail-s11.0-ms2000.prf",
#     "nail-s11.0-ms2147483647.prf",
# ]

# NAIL_S120 = [
#     "nail-s12.0-prog.prf",
#     "nail-s12.0-ms300.prf",
#     "nail-s12.0-ms2000.prf",
#     "nail-s12.0-ms2147483647.prf",
# ]

# MMSEQS_S57 = [
#     "mmseqs-s5.7-ms300.prf",
#     "mmseqs-s5.7-ms2000.prf",
#     "mmseqs-s5.7-ms2147483647.prf",
# ]

# MMSEQS_S75 = [
#     "mmseqs-s7.5-ms300.prf",
#     "mmseqs-s7.5-ms2000.prf",
#     "mmseqs-s7.5-ms2147483647.prf",
# ]

# MMSEQS_S90 = [
#     "mmseqs-s9.0-ms300.prf",
#     "mmseqs-s9.0-ms300.prf",
#     "mmseqs-s9.0-ms2000.prf",
# ]

# MMSEQS_S100 = [
#     "mmseqs-s10.0-ms300.prf",
#     "mmseqs-s10.0-ms2000.prf",
#     "mmseqs-s10.0-ms2147483647.prf",
# ]

# MMSEQS_S110 = [
#     "mmseqs-s11.0-ms300.prf",
#     "mmseqs-s11.0-ms2000.prf",
#     "mmseqs-s11.0-ms2147483647.prf",
# ]

# MMSEQS_S120 = [
#     "mmseqs-s12.0-ms300.prf",
#     "mmseqs-s12.0-ms2000.prf",
#     "mmseqs-s12.0-ms2147483647.prf",
# ]

NAIL_PROG = [
    # "nail-s5.7-prog.prf",
    # "nail-s7.5-prog.prf",
    "nail-s9.0-prog.prf",
    "nail-s10.0-prog.prf",
    "nail-s11.0-prog.prf",
    "nail-s12.0-prog.prf",
]

NAIL_300 = [
    # "nail-s5.7-ms300.prf",
    # "nail-s7.5-ms300.prf",
    "nail-s9.0-ms300.prf",
    "nail-s10.0-ms300.prf",
    "nail-s11.0-ms300.prf",
    "nail-s12.0-ms300.prf",
]

NAIL_2000 = [
    # "nail-s5.7-ms2000.prf",
    # "nail-s7.5-ms2000.prf",
    "nail-s9.0-ms2000.prf",
    "nail-s10.0-ms2000.prf",
    "nail-s11.0-ms2000.prf",
    "nail-s12.0-ms2000.prf",
]

NAIL_MAX = [
    # "nail-s5.7-ms2147483647.prf",
    # "nail-s7.5-ms2147483647.prf",
    "nail-s9.0-ms2147483647.prf",
    "nail-s10.0-ms2147483647.prf",
    "nail-s11.0-ms2147483647.prf",
    "nail-s12.0-ms2147483647.prf",
]

MMSEQS_300 = [
    # "mmseqs-s5.7-ms300.prf",
    "mmseqs-s7.5-ms300.prf",
    "mmseqs-s9.0-ms300.prf",
    "mmseqs-s10.0-ms300.prf",
    "mmseqs-s11.0-ms300.prf",
    "mmseqs-s12.0-ms300.prf",
]

MMSEQS_2000 = [
    # "mmseqs-s5.7-ms2000.prf",
    "mmseqs-s7.5-ms2000.prf",
    "mmseqs-s9.0-ms2000.prf",
    "mmseqs-s10.0-ms2000.prf",
    "mmseqs-s11.0-ms2000.prf",
    "mmseqs-s12.0-ms2000.prf",
]

MMSEQS_MAX = [
    # "mmseqs-s5.7-ms2147483647.prf",
    "mmseqs-s7.5-ms2147483647.prf",
    "mmseqs-s9.0-ms2147483647.prf",
    "mmseqs-s10.0-ms2147483647.prf",
    "mmseqs-s11.0-ms2147483647.prf",
    "mmseqs-s12.0-ms2147483647.prf",
]


def plot_line(pts, label, layout=None, label_layout=None):
    global ax

    pts.sort(key=lambda p: p.x)

    color = TOOL_COLORS[label.split()[0]]

    x = np.array([p.x for p in pts])
    y = np.array([p.y for p in pts])

    ax.plot(
        x, y,
        color=color,
        linestyle='--'
    )

    ax.scatter(
        x, y,
        color=color,
        linewidths=2,
        s=30,
    )

    if label_layout is None:
        o = (0, -20)
        ha = "center"
    else:
        o = label_layout[0]
        ha = label_layout[1]

    label_pt = (pts[0].x, pts[0].y)
    annotate(
        ax, label, label_pt, o, color,
        fontsize=12,
        rotation=0,
        linestyle='-',
        arrowstyle='-|>'
    )

    if layout is None:
        offsets = [(5, -5) for _ in range(len(pts))]
        ha_list = ["left" for _ in range(len(pts))]
    else:
        offsets = [lt[0] for lt in layout]
        ha_list = [lt[1] for lt in layout]

    for (p, o, ha) in zip(pts, offsets, ha_list):
        _, args = prefix_label(
            p.prefix,
            exclude=["ms", "prog"],
            reformat=[
                ("--mmseqs-s", "-s"),
                ("(sensitive)", ""),
                ("(default)", "")
            ]
        )
        args = " ".join(args)

        pt = (p.x, p.y)

        annotate(
            ax, f"{args}", pt, o, color,
            fontsize=8,
            rotation=0,
        )


if __name__ == "__main__":
    import sys
    filename = Path(sys.argv[1])

    points = []

    with open(filename) as f:
        lines = list(
            filter(
                lambda line: not line.startswith("#"),
                f.readlines()
            )
        )

        for line in lines:
            points.append(Point(line))

    n_300 = list(filter(lambda p: p.prefix in NAIL_300, points))
    n_2000 = list(filter(lambda p: p.prefix in NAIL_2000, points))
    n_max = list(filter(lambda p: p.prefix in NAIL_MAX, points))
    n_prog = list(filter(lambda p: p.prefix in NAIL_PROG, points))

    m_300 = list(filter(lambda p: p.prefix in MMSEQS_300, points))
    m_2000 = list(filter(lambda p: p.prefix in MMSEQS_2000, points))
    m_max = list(filter(lambda p: p.prefix in MMSEQS_MAX, points))

    global ax, fig
    fig, ax = axes()

    plot_line(
        n_prog, "nail (prog)",
        layout=[[(-5, 5), "left"] for _ in range(len(n_prog))],
        label_layout=[(0, -20), "center"]
    )

    plot_line(
        n_300, "nail (ms=300)",
        layout=[[(-5, 5), "left"] for _ in range(len(n_300))],
        label_layout=[(0, 20), "center"]
    )

    plot_line(
        n_2000, "nail (ms=2000)",
        layout=[[(5, -5), "left"] for _ in range(len(n_2000))],
        label_layout=[(0, -30), "center"]
    )

    plot_line(
        n_max, "nail (ms=max)",
        layout=[[(-5, 5), "left"] for _ in range(len(n_max))],
        label_layout=[(0, -15), "center"]
    )

    plot_line(
        m_300, "mmseqs (ms=300)",
        label_layout=[(0, 20), "center"]
    )

    plot_line(
        m_2000, "mmseqs (ms=2000)",
        label_layout=[(0, -20), "center"]
    )

    plot_line(
        m_max, "mmseqs (ms=max)",
        layout=[[(-5, 5), "left"] for _ in range(len(m_max))],
        label_layout=[(0, 20), "left"]
    )

    plt.title("")

    ax.set_xlabel("Runtime (seconds; log scale)")
    ax.set_xscale("log")
    ax.set_xlim(10.0, 4e3)

    ax.set_ylabel("Fraction of HMMER recall")
    ax.set_ylim(0.35, 1.025)

    ax.grid(True)

    plt.savefig(filename.with_suffix(".pdf"))
