#!/usr/bin/env python3
"""Where nail lost the hits hmmer found, as a flow left to right.

    plot.py stages.tbl --out figures/

hmmer's hits split into the ones nail found and the ones it missed. The
misses split by the stage nail lost them at, the ones scored under nail's
cutoff by whether hmmer's best domain clears its cutoff alone, the ones no
domain clears by how many domains hmmer reported, and the multi-domain ones
by the paper's classes for how their domains add up.

Every column fills the full height, split in proportion to its counts, and
the ribbons into it fan out from the one slice of the column before that it
splits. No slice is drawn thinner than one line of text, so the smallest are
over their share. Each slice is labelled inside the ribbon that runs into it.

The classes are read from cutoff-pairs.tbl beside stages.tbl.
"""

import argparse
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MPath

HMMER, NAIL = "#EE7733", "#EE3377"
BLUE, TEAL, RED, GREY = "#0077BB", "#009988", "#CC3311", "#BBBBBB"

STAGES = [
    ("lost_prefilter", "mmseqs prefilter"),
    ("lost_stopped", "prog seeding stopped"),
    ("lost_mmseqs", "mmseqs dropped"),
    ("lost_align", "nail filtered"),
    ("lost_cutoff", "under nail's cutoff"),
]

# in order of how far a class can be trusted
CLASSES = [
    ("strong_match", "strong"),
    ("colinear_match", "colinear"),
    ("pileup:full", "pileup: full"),
    ("pileup:mid", "pileup: mid"),
    ("pileup:fragment", "pileup: fragment"),
    ("shuffled_match", "shuffled"),
    ("weak_match", "weak"),
]


def style(scale=1.75):
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams.update({
        "font.size": mpl.rcParams["font.size"] * scale,
        "pdf.fonttype": 42,
        "axes.axisbelow": True,
        "grid.color": "0.9",
        "grid.linewidth": 0.8,
    })


def table(path):
    names, rows = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#="):
            continue
        if line.startswith("#"):
            f = line.lstrip("#").split()
            if f and not f[0].startswith("-"):
                names = f
            continue
        if line.split():
            rows.append(dict(zip(names, line.split())))
    return rows


def load(stages):
    """The one stages row as ints, plus derived numbers, and the pairs."""
    rows = table(stages)
    if len(rows) != 1:
        sys.exit(f"expected one stages row, got {len(rows)}")
    r = rows[0]
    d = {k: int(v) for k, v in r.items() if v.isdigit()}
    d["sens"] = float(r["sens"])
    d["missed"] = d["truth"] - d["kept"]
    assert d["missed"] == sum(d[k] for k, _ in STAGES)
    assert d["lost_cutoff"] == d["cutoff_domain"] + d["cutoff_nodomain"]
    assert d["cutoff_nodomain"] == d["nodomain_single"] + d["nodomain_multi"]
    pairs = table(Path(stages).with_name("cutoff-pairs.tbl"))
    known = {k for k, _ in CLASSES}
    multi = {}
    for p in pairs:
        if p["kind"] == "multi":
            if p["class"] not in known:
                sys.exit(f"unknown class {p['class']}")
            multi[p["class"]] = multi.get(p["class"], 0) + 1
    d["classes"] = [(name, multi[k]) for k, name in CLASSES if multi.get(k)]
    d["class_total"] = sum(n for _, n in d["classes"])
    if d["class_total"] != d["nodomain_multi"]:
        print(f"warning: classes sum to {d['class_total']:,}, stages.tbl says {d['nodomain_multi']:,}",
              file=sys.stderr)
    return d


def pct(n, of):
    p = 100 * n / of
    return f"{p:.1f}%" if p >= 0.1 else f"{p:.2f}%" if p >= 0.01 else "<0.01%"


def shades(hue, n):
    """n shades of one hue, the hue itself first and paler after."""
    rgb = np.array(mpl.colors.to_rgb(hue))
    return [tuple(rgb + (1 - rgb) * t) for t in np.linspace(0, 0.6, n)] if n > 1 else [hue]


def heights(counts, total, gap, least):
    """Each count's share of `total` less the gaps, no slice under `least`,
    the rest scaled down to make room."""
    fixed = set()
    while True:
        room = total - gap * (len(counts) - 1) - least * len(fixed)
        free = sum(n for i, n in enumerate(counts) if i not in fixed)
        hs = [least if i in fixed else room * n / free for i, n in enumerate(counts)]
        under = {i for i, h in enumerate(hs) if h < least} - fixed
        if not under:
            return hs
        fixed |= under


def draw(stages, out):
    style()
    d = load(stages)
    FS = mpl.rcParams["font.size"]

    W, H = 16.0, 8.0     # inches
    TOP = 0.7            # above the columns, for their headers
    NW = 0.22            # node width
    GAP = 0.12           # between the slices of a column
    LEAST = 0.24         # the thinnest slice, tall enough for one line

    fig = plt.figure(figsize=(W, H + TOP))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H + TOP)
    ax.axis("off")
    r = fig.canvas.get_renderer()

    # a column is its header, its slices top to bottom and which of them the
    # next column splits. the split slice is grey, and the rest are shades
    # of the column's own hue in slice order
    def column(header, hue, slices, cont):
        rest = iter(shades(hue, len(slices) - (cont is not None)))
        return header, [(name, n, GREY if i == cont else next(rest))
                        for i, (name, n) in enumerate(slices)], cont

    stage = [(name, d[key]) for key, name in STAGES]
    cols = [
        ("HMMER hits\n" + f"{d['truth']:,}", [("HMMER hits", d["truth"], HMMER)], 0),
        column("Found by nail", NAIL, [("found by nail", d["kept"]), ("missed", d["missed"])], 1),
        column("Stage lost at", BLUE, stage, len(stage) - 1),
        column("HMMER's best\ndomain", TEAL,
               [("a domain clears", d["cutoff_domain"]), ("no domain clears", d["cutoff_nodomain"])], 1),
        column("Domains HMMER\nreported", TEAL,
               [("1 domain", d["nodomain_single"]), ("2+ domains", d["nodomain_multi"])], 1),
        column("How 2+ domains\nadd up", RED, d["classes"], None),
    ]
    XS = np.linspace(0.3, W - 0.3 - NW, len(cols))

    def band(x0, y0a, y0b, x1, y1a, y1b, colour):
        m = (x0 + x1) / 2
        verts = [(x0, y0a), (m, y0a), (m, y1a), (x1, y1a), (x1, y1b), (m, y1b), (m, y0b), (x0, y0b), (x0, y0a)]
        codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.LINETO,
                 MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
        ax.add_patch(PathPatch(MPath(verts, codes), facecolor=colour, alpha=0.35, edgecolor="none"))

    def label(x, y0, y1, name, n, of):
        """Right-aligned against the node at x, two lines where they fit."""
        count = f"{n:,} ({pct(n, of)})"
        for text in (f"{name}\n{count}", f"{name}  {count}"):
            t = ax.text(x - 0.1, (y0 + y1) / 2, text, ha="right", va="center", fontsize=FS * 0.66,
                        color="0.1", linespacing=1.15)
            e = t.get_window_extent(r)
            if e.height <= (y1 - y0) * fig.dpi - 2:
                return
            t.remove()
        ax.text(x - 0.1, (y0 + y1) / 2, text, ha="right", va="center", fontsize=FS * 0.66, color="0.1")

    parent = None
    for x, (header, slices, cont) in zip(XS, cols):
        total = sum(n for _, n, _ in slices)
        hs = heights([n for _, n, _ in slices], H, GAP, LEAST)
        # the parent's slice, cut in proportion to the counts
        if parent is not None:
            px, py0, py1 = parent
            share = (py1 - py0) / total
        y = H
        top_down = []
        for (name, n, colour), h in zip(slices, hs):
            y0, y1 = y - h, y
            ax.add_patch(Rectangle((x, y0), NW, h, facecolor=colour, edgecolor="none"))
            if parent is not None:
                band(px + NW, py1 - share * n, py1, x, y0, y1, colour)
                py1 -= share * n
                label(x, y0, y1, name, n, total)
            top_down.append((y0, y1))
            y = y0 - GAP
        ax.text(x + NW / 2, H + 0.1, header, ha="center", va="bottom", fontsize=FS * 0.72, color="0.3",
                fontweight="bold")
        parent = None if cont is None else (x, *top_down[cont])

    path = out / "loss-decomp-hits.pdf"
    fig.savefig(path, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stages", help="the stages.tbl to plot")
    ap.add_argument("--out", default="figures", help="where the pdfs go")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    print(f"wrote {draw(args.stages, out)}")


if __name__ == "__main__":
    main()
