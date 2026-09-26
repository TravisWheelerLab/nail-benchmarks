"""Wall clock, and the percentage of perfect speedup reached, against thread
count, one column per unit.

Reads scaling.tbl as `thread-scaling parse` writes it: `#` lines are metadata
and the header, `-` is an empty cell. Writes thread-scaling.pdf into --out.
"""

import argparse
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt

# the first four slots of the validated categorical palette, in order, so an
# arm keeps its colour whichever of the others ran
ARMS = {
    "nail": "#2a78d6",
    "mmseqs": "#eb6834",
    "hmmer": "#1baf7a",
    "hmmer-split": "#eda100",
}
MARKERS = {"nail": "o", "mmseqs": "s", "hmmer": "^", "hmmer-split": "v"}


def read(path):
    names, rows = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            f = line.lstrip("#").split()
            if names is None and f and f[0] == "unit":
                names = f
            continue
        if not line.split():
            continue
        if names is None:
            raise SystemExit(f"no header in {path}")
        rows.append(dict(zip(names, line.split())))

    if not rows:
        raise SystemExit(f"no rows in {path}")
    return rows


def num(row, key):
    v = row.get(key, "-")
    return None if v == "-" else float(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("scaling")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = read(args.scaling)
    units = sorted({r["unit"] for r in rows})
    arms = [a for a in ARMS if any(r["arm"] == a for r in rows)]
    arms += sorted({r["arm"] for r in rows} - set(arms))

    fig, axes = plt.subplots(
        2, len(units), figsize=(4.2 * len(units), 7), squeeze=False, sharex=True
    )

    for col, unit in enumerate(units):
        top, bottom = axes[0][col], axes[1][col]
        threads_all = sorted({int(r["threads"]) for r in rows if r["unit"] == unit})

        for arm in arms:
            pts = sorted(
                (int(r["threads"]), num(r, "wall_s"), num(r, "ideal_s"), num(r, "pct_ideal"))
                for r in rows
                if r["unit"] == unit and r["arm"] == arm and num(r, "wall_s")
            )
            if not pts:
                continue
            t = [p[0] for p in pts]
            color = ARMS.get(arm, "#777777")
            kw = dict(
                color=color,
                marker=MARKERS.get(arm, "o"),
                markersize=6,
                linewidth=2,
                label=arm,
            )
            top.plot(t, [p[1] for p in pts], **kw)
            # perfect scaling from this arm's lowest rung, in its colour
            top.plot(t, [p[2] for p in pts], color=color, linestyle="--", linewidth=1)
            bottom.plot(t, [p[3] for p in pts], **kw)

        bottom.axhline(100.0, color="#999999", linestyle="--", linewidth=1)

        top.set_title(unit)
        top.set_xscale("log", base=2)
        top.set_yscale("log")
        bottom.set_xscale("log", base=2)
        bottom.set_ylim(0, 110)
        bottom.set_xticks(threads_all, [str(t) for t in threads_all])
        bottom.set_xlabel("threads")
        for ax in (top, bottom):
            ax.grid(True, color="#e5e5e5", linewidth=0.6)
            ax.set_axisbelow(True)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)

    axes[0][0].set_ylabel("wall clock (s); dashed = perfect scaling")
    axes[1][0].set_ylabel("% of perfect speedup")
    axes[0][0].legend(frameon=False, fontsize=9)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    path = out / "thread-scaling.pdf"
    fig.savefig(path)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
