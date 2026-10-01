#!/usr/bin/env python3
"""Where in nail's pipeline the ceiling's hits were lost.

    plot.py stages.tbl --out figures/

One panel per target corpus, one bar per arm. A bar is the hits the ceiling
kept that the arm did not, split by the stage that dropped them: never
seeded, seeded but absent from nail's table, or scored under the family's
cutoff. The bar is drawn to 100% of what was lost, so the arms compare by
where the loss sits rather than by how much there was, and how much there was
goes beside each bar with the arm's sensitivity.

The lost_ columns of stages.tbl are the stages, in the table's order, so a
column added there becomes a segment here; a column with no colour below is
an error rather than a hue picked at random.
"""

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

SCALE = 1.4
mpl.rcParams.update({"font.size": mpl.rcParams["font.size"] * SCALE})

# the palette the other benchmarks use, so figures from all of them sit together
TOL_BLUE = "#0077BB"
TOL_TEAL = "#009988"
TOL_RED = "#CC3311"

# a stage's column in stages.tbl, what to call it, and its colour.
# fixed per stage rather than cycled, so a figure with fewer stages
# keeps the same colour on each
STAGES = {
    "lost_seed": ("never seeded", TOL_BLUE),
    "lost_align": ("seeded, not aligned", TOL_TEAL),
    "lost_cutoff": ("under the family cutoff", TOL_RED),
}

# a segment narrower than this gets no label inside it: the legend
# carries its identity and the table its count
LABEL_FROM = 0.06


def table(path):
    """A `#`-headed table: its column names, its rows as dicts, and its `#=`
    lines as a key to the words after it, the header being the last `#` line
    whose first field is not a dash."""
    names, rows, meta = None, [], defaultdict(list)
    for line in Path(path).read_text().splitlines():
        if line.startswith("#="):
            f = line[2:].split()
            meta[f[0]].append(f[1:])
            continue
        if line.startswith("#"):
            f = line.lstrip("#").split()
            if f and not f[0].startswith("-"):
                names = f
            continue
        if line.split():
            rows.append(dict(zip(names, line.split())))
    return names, rows, meta


def draw(stages_path, out):
    names, rows, meta = table(stages_path)
    ceiling = meta["ceiling"][0][0]
    hits = {unit: int(n) for unit, n in meta["hits"]}
    stages = [n for n in names if n.startswith("lost_")]
    unknown = [n for n in stages if n not in STAGES]
    if unknown:
        raise SystemExit(f"no colour for {', '.join(unknown)}; add it to STAGES in {__file__}")

    units = list(dict.fromkeys(r["unit"] for r in rows))
    per_unit = {u: [r for r in rows if r["unit"] == u] for u in units}
    tallest = max(len(v) for v in per_unit.values())

    fig, axes = plt.subplots(
        len(units), 1, sharex=True, squeeze=False,
        figsize=(11, 1.4 + 0.85 * tallest * len(units)),
        layout="constrained",
    )
    for ax, unit in zip(axes[:, 0], units):
        arms = per_unit[unit]
        ys = range(len(arms))[::-1]
        for y, r in zip(ys, arms):
            lost = int(r["ceiling"]) - int(r["kept"])
            left = 0.0
            for stage in stages:
                n = int(r[stage])
                share = n / lost if lost else 0.0
                _, colour = STAGES[stage]
                # a surface-coloured edge is the gap between segments
                ax.barh(y, share, left=left, height=0.62, color=colour,
                        edgecolor="white", linewidth=1.5)
                if share >= LABEL_FROM:
                    ax.text(left + share / 2, y, f"{share:.0%}", ha="center",
                            va="center", color="white", fontsize=mpl.rcParams["font.size"] * 0.9)
                left += share
            # beside the bar: how much was lost, and what that leaves
            ax.text(1.03, y, f"{lost:,}", va="center", ha="left", color="0.2")
            ax.text(1.23, y, r["sens"], va="center", ha="left", color="0.2")

        ax.set_yticks(list(ys))
        ax.set_yticklabels([r["run"] for r in arms])
        ax.set_ylim(-0.6, len(arms) - 0.4)
        ax.set_title(f"{unit}    ceiling {ceiling}: {hits[unit]:,} hits",
                     loc="left", fontsize=mpl.rcParams["font.size"])
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", alpha=0.25)
        ax.set_axisbelow(True)

    top = axes[0, 0]
    top.text(1.03, len(per_unit[units[0]]) - 0.3, "hits lost", ha="left", va="bottom",
             color="0.4", fontsize=mpl.rcParams["font.size"] * 0.85)
    top.text(1.23, len(per_unit[units[0]]) - 0.3, "sensitivity", ha="left", va="bottom",
             color="0.4", fontsize=mpl.rcParams["font.size"] * 0.85)

    bottom = axes[-1, 0]
    bottom.set_xlim(0, 1.42)
    bottom.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    bottom.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
    bottom.set_xlabel("share of the ceiling's hits the arm lost, by stage")

    fig.legend(
        handles=[Patch(color=c, label=l) for l, c in (STAGES[s] for s in stages)],
        loc="outside lower center", ncol=len(stages), frameon=False,
    )
    fig.suptitle("loss-decomp: where nail loses the ceiling's hits", x=0.0, ha="left")

    path = out / "loss-decomp-stages.pdf"
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
