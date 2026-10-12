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
import re
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea, VPacker
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
from matplotlib.path import Path as MplPath

SCALE = 0.8
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
    "lost_prefilter": "mmseqs\nprefilter",
    "lost_stopped": "prog\nstopped",
    "lost_mmseqs": "mmseqs\ndropped",
    "lost_align": "nail\nfiltered",
    "lost_cutoff": "under\ncutoff",
}

# the paper's classes in cutoff-pairs.tbl's spelling and what to call each,
# in the paper's order
CLASSES = {
    "strong_match": "strong",
    "pileup:full": "pileup:\nfull",
    "pileup:mid": "pileup:\nmid",
    "pileup:fragment": "pileup:\nfragment",
    "colinear_match": "colinear",
    "shuffled_match": "shuffled",
    "weak_match": "weak",
    "-": "no positive\ndomain",
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


# the grid, in inches: one unit's four rows of two cells fill a letter
# page, each row a quarter of its height and each column half its width,
# less the margin round them, the gap between them and the title over them
PAGE_IN = (8.5, 11.0)
GAP_IN = 0.25
MARGIN_IN = 0.25
TITLE_IN = 0.4
CELL_W = (PAGE_IN[0] - 2 * MARGIN_IN - GAP_IN) / 2
CELL_H = (PAGE_IN[1] - 2 * MARGIN_IN - TITLE_IN - 3 * GAP_IN) / 4

# where a panel's axes go inside its cell, in inches from each edge: the
# title takes the top, and a panel with axes keeps room for its labels
FREE = (0.05, 0.4, 0.05, 0.05)
AXES = (0.55, 0.45, 0.1, 0.42)


def cell(fig, row, col):
    """The box of the cell at `row`, `col` of the grid, in figure
    coordinates, rows counted from the top under the title."""
    w, h = fig.get_size_inches()
    x0 = MARGIN_IN + col * (CELL_W + GAP_IN)
    y1 = h - MARGIN_IN - TITLE_IN - row * (CELL_H + GAP_IN)
    return mpl.transforms.Bbox.from_extents(x0 / w, (y1 - CELL_H) / h, (x0 + CELL_W) / w, y1 / h)


def axes_in(fig, box, margins):
    """Axes inside `box` with `margins` inches left, top, right and bottom
    of it free."""
    w, h = fig.get_size_inches()
    left, top, right, bottom = margins
    return fig.add_axes([
        box.x0 + left / w, box.y0 + bottom / h,
        box.width - (left + right) / w, box.height - (top + bottom) / h,
    ])


def frame(fig, box, title, colour=None):
    """A thin box around a panel, with its title inside at the top left.
    With `colour`, the count the title opens with is set in it: the colour
    of the arrow that brings those hits in from the panel above."""
    pad = 0.006
    box = box.expanded((box.width + 2 * pad) / box.width, (box.height + 2 * pad) / box.height)
    fig.add_artist(Rectangle(
        (box.x0, box.y0), box.width, box.height, transform=fig.transFigure, fill=False,
        edgecolor="0.7", linewidth=0.8, zorder=5,
    ))
    style = {"fontsize": mpl.rcParams["font.size"] * 0.9, "fontweight": "medium"}
    x, y = box.x0 + pad, box.y1 - pad
    first, rest = title.split("\n", 1) if "\n" in title else (title, None)
    count = re.search(r"\d[\d,]*", first) if colour else None
    if count is None:
        fig.text(x, y, title, ha="left", va="top", color=INK, linespacing=1.1, **style)
        return box

    # the first line in three pieces packed end to end, measured by
    # whichever renderer draws them, since the pdf's font metrics are not
    # the screen's
    pieces = [
        (first[:count.start()], INK),
        (count.group(), colour),
        (first[count.end():], INK),
    ]
    line = HPacker(children=[
        TextArea(text, textprops={"color": c, **style}) for text, c in pieces
    ], pad=0, sep=0)
    lines = [line] if rest is None else [line, TextArea(rest, textprops={"color": INK, **style})]
    fig.add_artist(AnchoredOffsetbox(
        loc="upper left", child=VPacker(children=lines, align="left", pad=0, sep=2),
        bbox_to_anchor=(x, y), bbox_transform=fig.transFigure, frameon=False, pad=0, borderpad=0,
    ))
    return box

    # the first line in three pieces, each starting where the last ended.
    # a piece's leading and trailing spaces are not drawn, so the gap at
    # each edge is stepped over by hand, a space being the width "a a"
    # has over "aa"
    renderer = fig.canvas.get_renderer()
    to_fig = fig.transFigure.inverted()

    def width(text):
        t = fig.text(0, 0, text, ha="left", **style)
        w = to_fig.transform_bbox(t.get_window_extent(renderer)).width
        t.remove()
        return w

    space = width("a a") - width("aa")
    pieces = [
        (first[:count.start()], INK),
        (count.group(), colour),
        (first[count.end():], INK),
    ]
    cx = x
    for text, c in pieces:
        cx += space * (len(text) - len(text.lstrip(" ")))
        t = fig.text(cx, y, text.strip(" "), ha="left", color=c, **style)
        extent = to_fig.transform_bbox(t.get_window_extent(renderer))
        cx = extent.x1 + space * (len(text) - len(text.rstrip(" ")))
    if rest:
        fig.text(x, y - extent.height * 1.1, rest, ha="left", color=INK, linespacing=1.1, **style)
    return box


# a chart panel's vertical layout, in its own units: the tallest bar is 1
# high over an axis at 0, and the names hang under the axis
CHART_BOTTOM = -0.22
CHART_TOP = 1.2


def chart(ax, bars, hue, size=0.7):
    """One bar per `(name, n)` to scale, in shades of `hue` from dark to
    light, each with its count over it and its name under it at `size` of
    the base font. Returns each bar's shade."""
    top = max(n for _, n in bars) or 1
    shades = ramp(hue, len(bars))
    for i, ((name, n), colour) in enumerate(zip(bars, shades)):
        ax.bar(i, n / top, width=0.7, color=colour, linewidth=0)
        ax.text(i, n / top + 0.02, f"{n:,}", ha="center", va="bottom", color="0.2",
                fontsize=mpl.rcParams["font.size"] * 0.85)

    ax.set_xlim(-0.5, len(bars) - 0.5)
    # no legend over a chart, so only the count over the tallest bar
    ax.set_ylim(CHART_BOTTOM, CHART_TOP)
    ax.set_xticks(range(len(bars)))
    ax.set_xticklabels([name for name, _ in bars], linespacing=1.0)
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_position(("data", 0))
    ax.spines["bottom"].set_color("0.6")
    # the size goes through tick_params: one set on the labels themselves
    # does not survive the ticks being rebuilt at draw
    ax.tick_params(axis="x", length=0, labelsize=mpl.rcParams["font.size"] * size)
    return shades



def pie(ax, segments, hue, angles):
    """One pie of slices to scale, each `(name, n)`, in shades of `hue`
    from dark to light, each with its share inside it where that fits and
    its name, count and share beside the rim. `angles` is where each
    slice's label goes, in degrees counterclockwise from the right, held
    inside the slice; the first slice is centred on its angle. Returns the
    labels, for the arrows that leave from them."""
    total = sum(n for _, n in segments)
    shades = ramp(hue, len(segments))
    start = angles[0] - 180 * segments[0][1] / total
    wedges, _ = ax.pie(
        [n for _, n in segments], colors=shades, startangle=start, radius=1.0,
        wedgeprops={"linewidth": 1.5, "edgecolor": "white"},
    )
    labels = []
    for (name, n), w, colour, want in zip(segments, wedges, shades, angles):
        share = n / total if total else 0.0
        mid = math.radians((w.theta1 + w.theta2) / 2)
        if share >= 0.08:
            ax.text(0.5 * math.cos(mid), 0.5 * math.sin(mid), f"{share:.0%}", ha="center",
                    va="center", color=ink_on(colour), fontsize=mpl.rcParams["font.size"] * 0.9)
        # the label's angle, taken round the circle to the slice's own
        # turn, stays a little inside the slice's edges
        margin = min(2.0, (w.theta2 - w.theta1) / 4)
        want = w.theta1 + (want - w.theta1) % 360
        a = math.radians(min(max(want, w.theta1 + margin), w.theta2 - margin))
        c, s_ = math.cos(a), math.sin(a)
        ax.plot([c, 1.2 * c], [s_, 1.2 * s_], color=ink(colour), linewidth=1.0,
                solid_capstyle="butt")
        ha = "left" if c > 0.3 else "right" if c < -0.3 else "center"
        va = "center" if abs(s_) < 0.3 else "bottom" if s_ > 0 else "top"
        labels.append(ax.text(
            1.28 * c, 1.28 * s_, f"{name}\n{n:,} ({share:.0%})", ha=ha, va=va, color=ink(colour),
            fontsize=mpl.rcParams["font.size"] * 0.8, linespacing=1.1,
        ))
    ax.set_xlim(-2.9, 2.9)
    ax.set_ylim(-1.25, 1.25)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    return labels



# the venn's vertical layout: the circles sit on the line at 0.5 with the
# larger one VENN_R high either side of it
VENN_R = 0.9
VENN_BOTTOM = -0.55
VENN_TOP = 1.55


def lens(r1, r2, d):
    """The area two circles of radii `r1` and `r2` share with their centres
    `d` apart."""
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2
    a = r1 * r1 * math.acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1))
    b = r2 * r2 * math.acos((d * d + r2 * r2 - r1 * r1) / (2 * d * r2))
    c = 0.5 * math.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    return a + b - c


def venn(ax, only_left, both, only_right, names, hue):
    """Two circles to scale, the left `only_left + both` in area and the
    right `both + only_right`, overlapping by `both`, in shades of `hue`
    from dark to light: left crescent, lens, right crescent. The lens is
    named inside it and each crescent beside it, with its count and share.
    Returns each region's `(left, right)` along the centre line and its
    shade, and the right crescent's label for an arrow to leave from."""
    total = only_left + both + only_right
    counts = [only_left, both, only_right]
    shades = ramp(hue, 3)
    scale = VENN_R / math.sqrt(max(only_left + both, both + only_right) / math.pi)
    r1 = scale * math.sqrt((only_left + both) / math.pi)
    r2 = scale * math.sqrt((both + only_right) / math.pi)

    # the distance between the centres that shares `both`, by bisection
    # on the lens area, which falls as the centres part
    want = scale * scale * both
    lo, hi = abs(r1 - r2), r1 + r2
    for _ in range(60):
        d = (lo + hi) / 2
        if lens(r1, r2, d) > want:
            lo = d
        else:
            hi = d
    d = (lo + hi) / 2
    y = 0.5
    c1, c2 = (-d / 2, y), (d / 2, y)

    ax.add_patch(Circle(c1, r1, color=shades[0], linewidth=0))
    ax.add_patch(Circle(c2, r2, color=shades[2], linewidth=0))
    overlap = ax.add_patch(Circle(c2, r2, color=shades[1], linewidth=0))
    overlap.set_clip_path(Circle(c1, r1, transform=ax.transData))
    # white edges over the fills are the gaps between regions
    for c, r in ((c1, r1), (c2, r2)):
        ax.add_patch(Circle(c, r, fill=False, edgecolor="white", linewidth=1.5))

    # each region's extent along the centre line
    spans = [
        (c1[0] - r1, c2[0] - r2),
        (c2[0] - r2, c1[0] + r1),
        (c1[0] + r1, c2[0] + r2),
    ]
    labels = [f"{name}\n{n:,}\n({n / total:.0%})" for name, n in zip(names, counts)]
    size = mpl.rcParams["font.size"] * 0.8
    # the lens is wide enough to name inside; the crescents are named
    # beside the circles, outside, in the crescent's ink
    (x0, x1), colour = spans[1], shades[1]
    ax.text((x0 + x1) / 2, y, labels[1].replace("\n(", " ("), ha="center", va="center",
            color=ink_on(colour), fontsize=size, linespacing=1.1)
    ax.text(spans[0][0] - 0.08, y, labels[0], ha="right", va="center",
                           color=ink(shades[0]), fontsize=size, linespacing=1.1)
    right = ax.text(spans[2][1] + 0.08, y, labels[2], ha="left", va="center",
                                   color=ink(shades[2]), fontsize=size, linespacing=1.1)

    ax.set_ylim(VENN_BOTTOM, VENN_TOP)
    ax.set_xlim(c1[0] - r1 - 0.2, c2[0] + r2 + 0.2)
    ax.axis("off")
    return [(span, shade) for span, shade in zip(spans, shades)], right


def round_off(fig, ax):
    """Widen a venn's x range to the aspect of the axes the layout gave it,
    so its circles come out round and the axes keep the whole cell. Called
    once the layout is frozen."""
    w, h = ax.get_window_extent(fig.canvas.get_renderer()).size
    y0, y1 = ax.get_ylim()
    x0, x1 = ax.get_xlim()
    xm, half = (x0 + x1) / 2, (y1 - y0) * (w / h) / 2
    ax.set_xlim(xm - half, xm + half)


def arrow(fig, tip, parent, box, colour):
    """An arrow from the label `tip` in one panel to the box of another,
    from the label's edge nearest the box: down, across and down again to
    the top edge of a box below, or straight sideways to a box in the same
    row, in figure coordinates. Placed once drawn, since a label's extent
    is only known then."""
    to_fig = fig.transFigure.inverted()
    b = to_fig.transform_bbox(tip.get_window_extent(fig.canvas.get_renderer()))
    xm, ym = (b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2
    if box.y1 <= b.y0:
        pad = box.width * 0.08
        xb = min(max(xm, box.x0 + pad), box.x1 - pad)
        # the parent's own box sits under the label, so the first leg
        # clears it before turning
        yc = (parent.y0 + box.y1) / 2
        points = [(xm, b.y0 - 0.003), (xm, yc), (xb, yc), (xb, box.y1)] if xb != xm \
            else [(xm, b.y0 - 0.003), (xb, box.y1)]
    elif (box.x0 + box.x1) / 2 < xm:
        points = [(b.x0 - 0.004, ym), (box.x1, ym)]
    else:
        points = [(b.x1 + 0.004, ym), (box.x0, ym)]
    fig.add_artist(FancyArrowPatch(
        path=MplPath(points), transform=fig.transFigure, arrowstyle="-|>", mutation_scale=14,
        color=colour, linewidth=1.0, shrinkA=0, shrinkB=0,
    ))


def histogram(ax, gaps, colour):
    """How far under the cutoff nail scored the pairs: quarter-bit bins up to
    one bit, whole bits past it, every bin drawn the same width."""
    if not gaps:
        ax.text(0.5, 0.5, "no pairs", ha="center", va="center", transform=ax.transAxes, color=colour)
        ax.axis("off")
        return
    edges = [0, 0.25, 0.5, 0.75] + list(range(1, max(2, math.ceil(max(gaps))) + 1))
    counts = [sum(1 for g in gaps if lo <= g < hi) for lo, hi in zip(edges, edges[1:])]
    slots = range(len(counts))
    ax.bar(slots, counts, width=1.0, align="edge", color=colour, edgecolor="white", linewidth=1.0)

    # a count over each bar reads only while the bars are wide; past a
    # dozen the axis carries the scale instead
    if len(counts) <= 12:
        for slot, n in zip(slots, counts):
            if n:
                ax.text(slot + 0.5, n, f"{n:,}", ha="center", va="bottom", color="0.2",
                        fontsize=mpl.rcParams["font.size"] * 0.85)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
    else:
        ax.tick_params(axis="y", labelsize=mpl.rcParams["font.size"] * 0.85)
        ax.spines["left"].set_color("0.6")

    gaps = sorted(gaps)
    median = gaps[len(gaps) // 2]
    ax.text(
        0.98, 0.95, f"{len(gaps):,} pairs, median {median:.1f} bits under",
        ha="right", va="top", transform=ax.transAxes, color=colour,
        fontsize=mpl.rcParams["font.size"] * 0.9,
    )

    # tick every edge while they fit, then the whole bits only, thinned
    ax.set_xlim(0, len(counts))
    labels = [f"{e:g}" for e in edges]
    if len(edges) <= 10:
        ticks = list(range(len(edges)))
    else:
        step = max(1, math.ceil((len(edges) - 4) / 8))
        ticks = [0, 2, 4] + list(range(4 + step, len(edges), step))
    ax.set_xticks(ticks)
    ax.set_xticklabels([labels[t] for t in ticks], fontsize=mpl.rcParams["font.size"] * 0.85)
    # the axis ends just over the tallest bin
    ax.set_ylim(0, max(counts) * 1.03)
    ax.set_xlabel("nail's score under the cutoff, bits")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def scatter(ax, pairs, colour):
    """hmmer's best domain score against its sequence score, with the line
    they would sit on were the sequence score the domain's."""
    if not pairs:
        ax.text(0.5, 0.5, "no pairs", ha="center", va="center", transform=ax.transAxes, color=colour)
        ax.axis("off")
        return
    xs = [float(p["domain"]) for p in pairs]
    ys = [float(p["hmmer"]) for p in pairs]
    lo, hi = min(xs + ys), max(xs + ys)
    pad = (hi - lo) * 0.05 or 1.0
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="0.7", linewidth=0.8)
    ax.text(
        hi + pad, hi + pad, "sequence = domain ", ha="right", va="top", color="0.5",
        fontsize=mpl.rcParams["font.size"] * 0.8, rotation=45, rotation_mode="anchor",
    )
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

    # four rows of two cells per unit: the venn beside the bar chart of the
    # stages; the histogram beside the pie of the cutoff loss; the pie of
    # the rest alone on the right; the scatter beside the bar chart of the
    # classes. Placed by hand, since a layout engine sizes a row to its
    # contents and these are to be equal
    rows_, cols_ = 4, 2 * len(units)
    fig = plt.figure(figsize=(2 * MARGIN_IN + cols_ * CELL_W + (cols_ - 1) * GAP_IN, PAGE_IN[1]))
    # each panel with its title and its cell
    panels = []
    arrows = []
    named = []
    venns = []

    for i, unit in enumerate(units):
        per_unit = [r for r in rows if r["unit"] == unit]
        if len(per_unit) != 1:
            raise SystemExit(f"{unit} has {len(per_unit)} nail runs in {stages_path}; one is drawn")
        r = per_unit[0]
        c = 2 * i
        cells = {(row, col): cell(fig, row, col) for row in range(rows_) for col in (c, c + 1)}
        first, second = axes_in(fig, cells[0, c], FREE), axes_in(fig, cells[0, c + 1], FREE)
        left, by_domain = axes_in(fig, cells[1, c], AXES), axes_in(fig, cells[1, c + 1], FREE)
        by_reported = axes_in(fig, cells[2, c + 1], FREE)
        points, classes = axes_in(fig, cells[3, c], AXES), axes_in(fig, cells[3, c + 1], FREE)

        truth, over, kept = int(r["truth"]), int(r["over"]), int(r["kept"])
        # hmmer's crescent on the right, facing the chart it leads to
        agreement, label = venn(first, over - kept, kept, truth - kept,
                                ("nail only", "both", "hmmer only"), HUES["agreement"])
        venns.append(first)
        panels.append((first, f"pairs over the family cutoff:\n{truth + over - kept:,} in either tool", cells[0, c]))

        lost = sum(int(r[s]) for s in stages)
        by_stage = chart(second, [(STAGES[s], int(r[s])) for s in stages], HUES["stages"])
        panels.append((second, f"the {lost:,} pairs only hmmer reported:\nby the stage nail lost them at", cells[0, c + 1]))

        under = int(r["lost_cutoff"])
        domain, nodomain = int(r["cutoff_domain"]), int(r["cutoff_nodomain"])
        # the slice that leads left is labelled on the left, the one that
        # leads down on the right
        domain_labels = pie(by_domain, [
            ("one domain\n≥ cutoff", domain),
            ("no domain\n≥ cutoff", nodomain),
        ], HUES["domain"], angles=(180, 0))
        panels.append((by_domain, f"the {under:,} under the cutoff:\nby hmmer's best domain", cells[1, c + 1]))

        pairs = [p for p in pair_rows if p["unit"] == unit and p["run"] == r["run"]]
        gaps = [float(p["cutoff"]) - float(p["nail"]) for p in pairs if p["kind"] == "domain"]
        histogram(left, gaps, ramp(HUES["gap"], 1)[0])
        panels.append((left, f"the {domain:,} one domain clears:\nnail's score under the cutoff", cells[1, c]))

        single, multi = int(r["nodomain_single"]), int(r["nodomain_multi"])
        reported_labels = pie(by_reported, [
            ("one reported\ndomain", single),
            ("2+ reported\ndomains", multi),
        ], HUES["reported"], angles=(180, 0))
        panels.append((by_reported, f"the {nodomain:,} no domain clears:\nby how many domains hmmer reported", cells[2, c + 1]))

        scatter(points, [p for p in pairs if p["kind"] == "single"], ramp(HUES["scatter"], 1)[0])
        panels.append((points, f"the {single:,} with one domain:\nits score against the sequence's", cells[3, c]))

        counts = {}
        for p in pairs:
            if p["kind"] == "multi":
                counts[p["class"]] = counts.get(p["class"], 0) + 1
        unknown = [k for k in counts if k not in CLASSES]
        if unknown:
            raise SystemExit(f"no colour for class {', '.join(unknown)}; add it to CLASSES in {__file__}")
        chart(classes, [
            (CLASSES[k], counts.get(k, 0)) for k in CLASSES if counts.get(k, 0)
        ], HUES["classes"], size=0.6)
        panels.append((classes, f"the {multi:,} with 2+ domains:\nby how they add up to hmmer's score", cells[3, c + 1]))

        arrows.append((first, label, second, ink(agreement[2][1])))
        # the arrow leaves from the under-cutoff bar's name, which is only
        # there once drawn
        shade = by_stage[stages.index("lost_cutoff")]
        named.append((second, stages.index("lost_cutoff"), by_domain, ink(shade)))
        shades = ramp(HUES["domain"], 2)
        arrows.append((by_domain, domain_labels[0], left, ink(shades[0])))
        arrows.append((by_domain, domain_labels[1], by_reported, ink(shades[1])))
        shades = ramp(HUES["reported"], 2)
        arrows.append((by_reported, reported_labels[0], points, ink(shades[0])))
        arrows.append((by_reported, reported_labels[1], classes, ink(shades[1])))

        if len(units) > 1:
            first.set_title(unit, loc="left", fontsize=mpl.rcParams["font.size"])

    w, h = fig.get_size_inches()
    fig.suptitle("loss-decomp: pairs over the cutoff, and where nail lost hmmer's",
                 x=MARGIN_IN / w, y=1 - MARGIN_IN / h, ha="left", va="top")

    # the arrows leave from labels whose extent is only known once drawn
    fig.canvas.draw()
    for ax in venns:
        round_off(fig, ax)
    for ax, i, below, colour in named:
        arrows.append((ax, ax.get_xticklabels()[i], below, colour))
    boxes = {ax: box for ax, _, box in panels}
    # a child's title opens with the count its arrow brings in, in the
    # arrow's colour
    brought = {below: colour for _, _, below, colour in arrows}
    for ax, title, box in panels:
        frame(fig, box, title, brought.get(ax))
    for above, tip, below, colour in arrows:
        arrow(fig, tip, boxes[above], boxes[below], colour)

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
