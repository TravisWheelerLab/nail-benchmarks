#!/usr/bin/env python3
"""What each seeding arm found, against what it cost.

    plot.py stages.tbl ledger.tbl --out figures/

One panel per target corpus. The two seeding policies sweep different knobs --
static moves --mmseqs-max-seqs alone, prog moves --prog-n against --prog-f --
so they share no axis of their own. What they do share is the bill: an arm is
a seeding plus the alignment that replays it, and sensitivity per second is the
question a reader actually has.

The y axis is what an arm missed rather than what it found, on a log scale.
Sensitivity 0.9103 against 0.9093 is invisible on an axis that also has to hold
--mmseqs-max-seqs 200 at 0.6246; the same pair as 8.97% missed against 9.07% is
legible, and the arm that misses 37.5% sits four times up the axis.

The stage decomposition is not drawn. Loss at the cloud and alignment stages
runs to tens of pairs against tens of thousands lost at seeding, so a stacked
bar would be one visible colour and a sentence says it better.
"""

import argparse
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt

SCALE = 1.5
mpl.rcParams.update({"font.size": mpl.rcParams["font.size"] * SCALE})

# the palette the other benchmarks use, so figures from all of them sit together
TOL_RED = "#CC3311"
TOL_TEAL = "#009988"


def rows(path):
    """A `#`-headed table as dicts, the header being the last `#` line whose
    first field is not a dash."""
    names, out = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            f = line.lstrip("#").split()
            if f and not f[0].startswith("-"):
                names = f
            continue
        if line.split():
            out.append(dict(zip(names, line.split())))
    return out


def costs(ledger):
    """Seconds per (unit, arm): the alignment, plus the seeding it replayed.

    A run's ledger row already sums the shards it covered, but these pipelines
    write one row per unit, so the pair is what keys it.
    """
    out = {}
    for r in rows(ledger):
        wall = float(r["wall(s)"]) if r["wall(s)"] != "-" else 0.0
        arm = r["seeds"] if r["seeds"] != "-" else r["name"]
        if arm == "-" or r["tool"] == "hmmer":
            continue
        out[(r["shard"], arm)] = out.get((r["shard"], arm), 0.0) + wall
    return out


def pareto(points):
    """The arms nothing else beats on both cost and sensitivity."""
    front, best = [], -1.0
    for p in sorted(points, key=lambda p: (p[0], -p[1])):
        if p[1] > best:
            front.append(p)
            best = p[1]
    return front


def draw(stages, ledger, out):
    cost = costs(ledger)
    at = {}
    for r in rows(stages):
        at.setdefault(r["unit"], []).append((r["run"], float(r["sens"])))

    units = sorted(at)
    fig, axes = plt.subplots(1, len(units), figsize=(9 * len(units), 8),
                             constrained_layout=True)
    axes = [axes] if len(units) == 1 else list(axes)

    for ax, unit in zip(axes, units):
        pts = []
        for arm, sens in at[unit]:
            seconds = cost.get((unit, arm))
            if seconds is None:
                continue
            pts.append((seconds, sens, arm))

        # static sweeps one knob and joins into a line; prog sweeps two and is
        # a cloud, so the shapes say which policy a point belongs to
        line = sorted(p for p in pts if p[2].startswith("static"))
        cloud = [p for p in pts if not p[2].startswith("static")]

        miss = lambda s: max(1.0 - s, 1e-4)

        ax.plot([p[0] for p in line], [miss(p[1]) for p in line], "-o",
                color=TOL_RED, markersize=11, linewidth=2.5, label="static", zorder=3)
        ax.scatter([p[0] for p in cloud], [miss(p[1]) for p in cloud],
                   s=110, color=TOL_TEAL, label="prog", zorder=3)

        front = pareto([(p[0], p[1]) for p in pts])
        ax.plot([p[0] for p in front], [miss(p[1]) for p in front], "--",
                color="0.4", linewidth=1.5, zorder=2, label="pareto front")

        # only the front is labelled: the arms behind it are the ones a reader
        # has no decision to make about, and twelve labels in one cluster is
        # nothing anyone can read
        on_front = {(round(c, 3), round(s, 6)) for c, s in front}
        for seconds, sens, arm in pts:
            if (round(seconds, 3), round(sens, 6)) not in on_front:
                continue
            ax.annotate(arm.replace("static-ms", "ms").replace("prog-", ""),
                        (seconds, miss(sens)), textcoords="offset points",
                        xytext=(8, 4), fontsize=10, color="0.2")

        ax.set_yscale("log")
        ax.set_title(unit)
        ax.set_xlabel("wall clock for the arm (s)   seeding + alignment")
        ax.set_ylabel("missed, of what hmmer found")
        ax.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda v, _: f"{v * 100:g}%")
        )
        ax.grid(alpha=0.25, which="both")
        ax.legend(loc="upper right")

    fig.suptitle("loss-decomp: what each seeding arm missed, against what it cost")

    path = out / "seeding.pdf"
    fig.savefig(path)
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stages", help="the stages.tbl to plot")
    ap.add_argument("ledger", help="the ledger.tbl its costs come from")
    ap.add_argument("--out", default="figures", help="where the pdfs go")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    print(f"wrote {draw(args.stages, args.ledger, out)}")


if __name__ == "__main__":
    main()
