#!/usr/bin/env python3
"""The pairs hmmer and nail report over the family cutoff, where nail lost
the ones only hmmer did, and what the ones it scored under the cutoff look
like in hmmer's domains.

    plot.py stages.tbl --out figures/

Per unit, a tree of boxed panels. A segment of one panel that the next
panel splits further is marked with a brace under it in the segment's
colour, its hit count under the brace, and an arrow in the same colour
runs from under the count to the panel that splits it, bending at right
angles.

The first panel is every pair either tool reported at or above the family's
cutoff: the pairs only hmmer reported, the pairs both did, the pairs only
nail did. The second is the pairs only hmmer reported, by the stage nail
lost them at: not in the mmseqs prefilter list, in it but past the rank
where prog seeding stopped for the query so never aligned, aligned by mmseqs
and dropped, seeded but dropped by nail's cloud search or forward filter, or
scored under the family's cutoff. The third is the pairs scored under the
cutoff, by whether one of hmmer's domains clears the cutoff on its own.
Under that, two panels: the pairs a domain clears open out into a histogram
of how far under the cutoff nail scored them, and the pairs no domain clears
into a bar of how many domains hmmer reported for them, one or more. Under
that, two more: the one-domain pairs into a scatter of hmmer's domain score
against its sequence score, which is the score hmmer's forward pass gave
the sequence over what it attributed to the domain, and the multi-domain
pairs into a bar of the paper's classes for how their domains add up. The
per-pair numbers and classes are read from cutoff-pairs.tbl beside
stages.tbl.

Every panel carries its own title and one hue: a bar's segments are that
hue from dark to light in bar order, and the brace, count and arrow under
a segment take its shade. Every segment is to scale within its bar,
carries its share inside it where that fits, and is named with its count
and share in the panel's legend over the bar. The histogram and the
scatter have the only axes.

The lost_ columns of stages.tbl are the stages, in the table's order, so a
column added there becomes a segment here; a column with no colour below is
an error rather than a hue picked at random.
"""

import argparse
import colorsys
import math
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Patch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath

SCALE = 1.4
mpl.rcParams.update({"font.size": mpl.rcParams["font.size"] * SCALE})

# one hue per panel, from the palette the other benchmarks use, so figures
# from all of them sit together
HUES = {
    "agreement": "#0077BB",
    "stages": "#CC3311",
    "domain": "#AA4499",
    "gap": "#AA4499",
    "reported": "#117A65",
    "scatter": "#117A65",
    "classes": "#EE7733",
}

# a stage's column in stages.tbl and what to call it, in pipeline order
STAGES = {
    "lost_prefilter": "pre-filtered (mmseqs)",
    "lost_stopped": "not aligned (prog stopped)",
    "lost_mmseqs": "aligned, dropped (mmseqs)",
    "lost_align": "filtered (nail)",
    "lost_cutoff": "under cutoff",
}

# the paper's classes in cutoff-pairs.tbl's spelling and what to call each,
# in the paper's order
CLASSES = {
    "strong_match": "strong",
    "pileup:full": "pileup:full",
    "pileup:mid": "pileup:mid",
    "pileup:fragment": "pileup:fragment",
    "colinear_match": "colinear",
    "shuffled_match": "shuffled",
    "weak_match": "weak",
    "-": "no positive domain",
}


def ramp(hue, n):
    """`n` shades of one hue from dark to light, as hex."""
    h, _, s = colorsys.rgb_to_hls(*mpl.colors.to_rgb(hue))
    s = min(s, 0.75)
    if n == 1:
        return [mpl.colors.to_hex(colorsys.hls_to_rgb(h, 0.42, s))]
    return [
        mpl.colors.to_hex(colorsys.hls_to_rgb(h, 0.3 + 0.42 * i / (n - 1), s))
        for i in range(n)
    ]


def ink(shade):
    """A shade dark enough to read as a line or text: the same hue, no
    lighter than the middle of the ramp."""
    h, l, s = colorsys.rgb_to_hls(*mpl.colors.to_rgb(shade))
    return mpl.colors.to_hex(colorsys.hls_to_rgb(h, min(l, 0.4), min(s, 0.6)))

# the width of one unit's panel in inches
WIDTH_IN = 9.0

INK = "0.35"


def ink_on(colour):
    """White or near-black, whichever clears the fill."""
    r, g, b = (int(colour[i : i + 2], 16) / 255 for i in (1, 3, 5))
    return "white" if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.5 else "0.15"


def table(path):
    """A `#`-headed table's rows as dicts, the header being the last `#` line
    whose first field is not a dash."""
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
    return names, rows


def outside_layout(artist):
    """The artist as drawn but not laid out: the layout engine sizes the axes
    to their contents, and a label that may hang past the axes' edge must
    not shrink them."""
    artist.set_in_layout(False)
    return artist


def extents(fig, panels):
    """Each panel's box in figure coordinates, once the layout is frozen: its
    whole extent, tick and axis labels included, widened to the row's
    tallest and the column's widest so the boxes line up as a grid.
    `panels` are `(ax, row, column)`, the column `None` for a full-width
    panel."""
    to_fig = fig.transFigure.inverted()
    renderer = fig.canvas.get_renderer()
    tight = {ax: to_fig.transform_bbox(ax.get_tightbbox(renderer)) for ax, _, _ in panels}

    rows, cols = {}, {}
    for ax, row, col in panels:
        b = tight[ax]
        y0, y1 = rows.get(row, (b.y0, b.y1))
        rows[row] = (min(y0, b.y0), max(y1, b.y1))
        x0, x1 = cols.get(col, (b.x0, b.x1))
        cols[col] = (min(x0, b.x0), max(x1, b.x1))
    # a full-width panel spans every column
    full = (min(x0 for x0, _ in cols.values()), max(x1 for _, x1 in cols.values()))

    return {
        ax: mpl.transforms.Bbox.from_extents(
            cols[col][0] if col is not None else full[0], rows[row][0],
            cols[col][1] if col is not None else full[1], rows[row][1],
        )
        for ax, row, col in panels
    }


def frame(fig, box, title):
    """A thin box around a panel, with its title inside at the top left."""
    pad = 0.006
    box = box.expanded((box.width + 2 * pad) / box.width, (box.height + 2 * pad) / box.height)
    fig.add_artist(Rectangle(
        (box.x0, box.y0), box.width, box.height, transform=fig.transFigure, fill=False,
        edgecolor="0.7", linewidth=0.8, zorder=5,
    ))
    fig.text(box.x0 + pad, box.y1 - pad, title, ha="left", va="top", color=INK,
             fontsize=mpl.rcParams["font.size"] * 0.9, fontweight="medium", linespacing=1.1)
    return box


def headroom(ax, lines):
    """Blank title lines over a panel with axes, so its box has room for the
    title that frame() writes at the box's top left."""
    ax.set_title("\n".join([" "] * lines), loc="left", fontsize=mpl.rcParams["font.size"] * 0.9)


# a bar panel's vertical layout, in the bar's own units: the bar is 0.6
# thick at 0, a brace and its count go under it to -1.0, and the legend
# goes over it
BAR_BOTTOM = -1.15
BAR_TOP = 1.9


def bar(ax, segments, hue, columns=3, stretch=1.0):
    """One bar of segments to scale, each `(name, n)`, in shades of `hue`
    from dark to light, each with its share inside it where that fits, and
    every one named with its count and share in a legend over the bar,
    `columns` wide. `stretch` is how much taller than a standard bar panel
    the axes are, so the bar comes out the same thickness. Returns each
    segment's `(left, right)` in the bar's own units and its shade."""
    total = sum(n for _, n in segments)
    spans = []
    left = 0
    handles = []
    for (name, n), colour in zip(segments, ramp(hue, len(segments))):
        share = n / total if total else 0.0
        # a surface-coloured edge is the gap between segments
        ax.barh(0, n, left=left, height=0.6, color=colour, edgecolor="white", linewidth=1.5)
        if share >= 0.08:
            ax.text(left + n / 2, 0, f"{share:.0%}", ha="center", va="center", color=ink_on(colour),
                    fontsize=mpl.rcParams["font.size"] * 0.9)
        handles.append(Patch(color=colour, label=f"{name}  {n:,}  ({share:.0%})"))
        spans.append(((left, left + n), colour))
        left += n

    # a taller panel's extra room all goes over the bar, where the legend is
    extra = (BAR_TOP - BAR_BOTTOM) * (stretch - 1)
    ax.set_xlim(-total * 0.01, total * 1.01)
    ax.set_ylim(BAR_BOTTOM, BAR_TOP + extra)
    ax.axis("off")
    # anchored just over the bar, in the bar's own units
    top = (0.4 - BAR_BOTTOM) / ((BAR_TOP + extra) - BAR_BOTTOM)
    outside_layout(ax.legend(
        handles=handles, loc="lower left", bbox_to_anchor=(0, top), ncol=columns,
        frameon=False, handlelength=1.0, handleheight=1.0, columnspacing=1.4,
        borderaxespad=0, fontsize=mpl.rcParams["font.size"] * 0.85,
    ))
    return spans


def brace(ax, span, count, colour):
    """A curly brace under one segment of a bar in the segment's colour, its
    tip at the middle and the segment's hit count under the tip. Returns
    where an arrow leaves from under the count, in data coordinates."""
    x0, x1 = span
    total = ax.get_xlim()[1] / 1.01
    y, h = -0.36, 0.14
    w = min((x1 - x0) / 4, total * 0.012)
    xm = (x0 + x1) / 2
    k = 0.55
    verts = [
        (x0, y),
        (x0, y - k * h), (x0 + (1 - k) * w, y - h), (x0 + w, y - h),
        (xm - w, y - h),
        (xm - (1 - k) * w, y - h), (xm, y - h - k * h), (xm, y - 2 * h),
        (xm, y - h - k * h), (xm + (1 - k) * w, y - h), (xm + w, y - h),
        (x1 - w, y - h),
        (x1 - (1 - k) * w, y - h), (x1, y - k * h), (x1, y),
    ]
    codes = [MplPath.MOVETO] + [MplPath.CURVE4] * 3 + [MplPath.LINETO] + [MplPath.CURVE4] * 6 \
        + [MplPath.LINETO] + [MplPath.CURVE4] * 3
    outside_layout(ax.add_patch(PathPatch(
        MplPath(verts, codes), fill=False, edgecolor=colour, linewidth=1.0, clip_on=False,
    )))
    outside_layout(ax.text(
        xm, y - 2 * h - 0.06, f"{count:,} hits", ha="center", va="top", color=colour,
        fontsize=mpl.rcParams["font.size"] * 0.85,
    ))
    return (xm, y - 2 * h - 0.42)


def arrow(fig, above, tip, parent, box, colour):
    """An arrow from under a brace's count in the panel above to the top
    edge of the box below, down, across and down again, in figure
    coordinates. Placed once the layout is frozen, since the boxes are only
    known then."""
    to_fig = fig.transFigure.inverted()
    (xa, ya) = to_fig.transform(above.transData.transform(tip))
    pad = box.width * 0.08
    xb = min(max(xa, box.x0 + pad), box.x1 - pad)
    yb = box.y1
    # the parent's own box sits under the brace, so the first leg clears it
    # before turning
    yc = (parent.y0 + yb) / 2
    points = [(xa, ya), (xa, yc), (xb, yc), (xb, yb)] if xb != xa else [(xa, ya), (xb, yb)]
    fig.add_artist(FancyArrowPatch(
        path=MplPath(points), transform=fig.transFigure, arrowstyle="-|>", mutation_scale=14,
        color=colour, linewidth=1.0, shrinkA=0, shrinkB=0,
    ))


def histogram(ax, gaps, colour):
    """How far under the cutoff nail scored the pairs, in half-bit bins, each
    bin counted above it."""
    headroom(ax, 2)
    if not gaps:
        ax.text(0.5, 0.5, "no pairs", ha="center", va="center", transform=ax.transAxes, color=colour)
        ax.axis("off")
        return
    width = 0.5
    edges = [i * width for i in range(max(1, math.ceil(max(gaps) / width)) + 1)]
    counts = [sum(1 for g in gaps if lo <= g < lo + width) for lo in edges[:-1]]
    ax.bar(edges[:-1], counts, width=width, align="edge", color=colour, edgecolor="white", linewidth=1.0)
    for lo, n in zip(edges[:-1], counts):
        if n:
            ax.text(lo + width / 2, n, f"{n:,}", ha="center", va="bottom", color="0.2",
                    fontsize=mpl.rcParams["font.size"] * 0.85)

    gaps = sorted(gaps)
    median = gaps[len(gaps) // 2]
    outside_layout(ax.text(
        0.98, 0.95, f"{len(gaps):,} pairs, median {median:.1f} bits under",
        ha="right", va="top", transform=ax.transAxes, color=colour,
        fontsize=mpl.rcParams["font.size"] * 0.9,
    ))

    ax.set_xlim(0, edges[-1])
    ax.set_ylim(0, max(counts) * 1.4)
    ax.set_xticks(edges[::2] if len(edges) > 8 else edges)
    ax.set_xlabel("nail's score under the cutoff, bits")
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)


def scatter(ax, pairs, colour):
    """hmmer's best domain score against its sequence score, with the line
    they would sit on were the sequence score the domain's."""
    headroom(ax, 2)
    if not pairs:
        ax.text(0.5, 0.5, "no pairs", ha="center", va="center", transform=ax.transAxes, color=colour)
        ax.axis("off")
        return
    xs = [float(p["domain"]) for p in pairs]
    ys = [float(p["hmmer"]) for p in pairs]
    lo, hi = min(xs + ys), max(xs + ys)
    pad = (hi - lo) * 0.05 or 1.0
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="0.7", linewidth=0.8)
    outside_layout(ax.text(
        hi + pad, hi + pad, "sequence = domain ", ha="right", va="top", color="0.5",
        fontsize=mpl.rcParams["font.size"] * 0.8, rotation=45, rotation_mode="anchor",
    ))
    ax.scatter(xs, ys, s=12, color=colour, alpha=0.6, linewidths=0)
    ax.text(0.03, 0.95, f"{len(pairs):,} pairs", ha="left", va="top", transform=ax.transAxes,
            color=colour, fontsize=mpl.rcParams["font.size"] * 0.9)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad * 2)
    ax.set_xlabel("hmmer's domain score, bits")
    ax.set_ylabel("hmmer's sequence score, bits")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def draw(stages_path, out):
    names, rows = table(stages_path)
    stages = [n for n in names if n.startswith("lost_")]
    unknown = [n for n in stages if n not in STAGES]
    if unknown:
        raise SystemExit(f"no colour for {', '.join(unknown)}; add it to STAGES in {__file__}")
    units = list(dict.fromkeys(r["unit"] for r in rows))

    _, pair_rows = table(Path(stages_path).with_name("cutoff-pairs.tbl"))

    # five rows per unit: three bars the width of the unit, then two panels
    # side by side under the halves of the third bar, then two more under
    # the halves of the right one of those
    # the bottom two rows hold a panel with axes beside a bar, and are as
    # tall as the axes need; the bar beside each is stretched to match
    heights = [1, 1.35, 1, 2.0, 2.4]
    fig = plt.figure(figsize=(WIDTH_IN * len(units), 15), layout="constrained")
    fig.get_layout_engine().set(hspace=0.22, wspace=0.1, rect=(0.012, 0.008, 0.976, 0.975))
    gs = fig.add_gridspec(5, 2 * len(units), height_ratios=heights)
    # each panel with its title, row and column, the column None for a
    # full-width one
    panels = []
    arrows = []

    for i, unit in enumerate(units):
        per_unit = [r for r in rows if r["unit"] == unit]
        if len(per_unit) != 1:
            raise SystemExit(f"{unit} has {len(per_unit)} nail runs in {stages_path}; one is drawn")
        r = per_unit[0]
        c = 2 * i
        first, second, third = (fig.add_subplot(gs[row, c:c + 2]) for row in range(3))
        left, right = fig.add_subplot(gs[3, c]), fig.add_subplot(gs[3, c + 1])
        points, classes = fig.add_subplot(gs[4, c]), fig.add_subplot(gs[4, c + 1])

        truth, over, kept = int(r["truth"]), int(r["over"]), int(r["kept"])
        agreement = bar(first, [
            ("hmmer only", truth - kept),
            ("both", kept),
            ("nail only", over - kept),
        ], HUES["agreement"])
        panels.append((first, f"pairs over the family cutoff, {truth + over - kept:,} in either tool", 0, None))

        lost = sum(int(r[s]) for s in stages)
        by_stage = bar(second, [(STAGES[s], int(r[s])) for s in stages], HUES["stages"],
                       columns=2, stretch=heights[1] / heights[0])
        panels.append((second, f"the {lost:,} pairs only hmmer reported, by the stage nail lost them at", 1, None))

        under = int(r["lost_cutoff"])
        domain, nodomain = int(r["cutoff_domain"]), int(r["cutoff_nodomain"])
        by_domain = bar(third, [
            ("one domain ≥ cutoff", domain),
            ("no domain ≥ cutoff", nodomain),
        ], HUES["domain"])
        panels.append((third, f"the {under:,} pairs nail scored under the cutoff, by hmmer's best domain", 2, None))

        pairs = [p for p in pair_rows if p["unit"] == unit and p["run"] == r["run"]]
        gaps = [float(p["cutoff"]) - float(p["nail"]) for p in pairs if p["kind"] == "domain"]
        histogram(left, gaps, ramp(HUES["gap"], 1)[0])
        panels.append((left, f"the {domain:,} one domain clears:\nnail's score under the cutoff", 3, c))

        single, multi = int(r["nodomain_single"]), int(r["nodomain_multi"])
        by_reported = bar(right, [
            ("one reported domain", single),
            ("2+ reported domains", multi),
        ], HUES["reported"], columns=1, stretch=heights[3] / heights[0])
        panels.append((right, f"the {nodomain:,} no domain clears:\nby how many domains hmmer reported", 3, c + 1))

        scatter(points, [p for p in pairs if p["kind"] == "single"], ramp(HUES["scatter"], 1)[0])
        panels.append((points, f"the {single:,} with one domain:\nits score against the sequence's", 4, c))

        counts = {}
        for p in pairs:
            if p["kind"] == "multi":
                counts[p["class"]] = counts.get(p["class"], 0) + 1
        unknown = [k for k in counts if k not in CLASSES]
        if unknown:
            raise SystemExit(f"no colour for class {', '.join(unknown)}; add it to CLASSES in {__file__}")
        bar(classes, [
            (CLASSES[k], counts.get(k, 0)) for k in CLASSES if counts.get(k, 0)
        ], HUES["classes"], columns=1, stretch=heights[4] / heights[0])
        panels.append((classes, f"the {multi:,} with 2+ domains:\nby how they add up to hmmer's score", 4, c + 1))

        for above, (span, shade), count, below in (
            (first, agreement[0], lost, second),
            (second, by_stage[stages.index("lost_cutoff")], under, third),
            (third, by_domain[0], domain, left),
            (third, by_domain[1], nodomain, right),
            (right, by_reported[0], single, points),
            (right, by_reported[1], multi, classes),
        ):
            arrows.append((above, brace(above, span, count, ink(shade)), below, ink(shade)))

        if len(units) > 1:
            first.set_title(unit, loc="left", fontsize=mpl.rcParams["font.size"])

    fig.suptitle("loss-decomp: pairs over the cutoff, and where nail lost hmmer's",
                 x=0.0, ha="left")

    # the boxes and arrows are figure-level shapes placed off the final
    # positions of the panels, so the layout is run once and then frozen
    # before they go on
    fig.canvas.draw()
    fig.set_layout_engine("none")
    boxes = extents(fig, [(ax, row, col) for ax, _, row, col in panels])
    for ax, title, _, _ in panels:
        frame(fig, boxes[ax], title)
    for above, tip, below, colour in arrows:
        arrow(fig, above, tip, boxes[above], boxes[below], colour)

    path = out / "loss-decomp-hits.pdf"
    fig.savefig(path)
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
