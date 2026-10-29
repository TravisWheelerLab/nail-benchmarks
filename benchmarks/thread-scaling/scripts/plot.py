"""Wall clock, and the percentage of perfect speedup reached, against thread
count: one figure each, with a unit's arms drawn together.

Reads scaling.tbl as `thread-scaling parse` writes it: `#` lines are metadata
and the header, `-` is an empty cell. Writes thread-scaling-wall.pdf and
thread-scaling-efficiency.pdf into --out, with the unit in the name when the
table holds more than one.
"""

import argparse
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt

# pid's canvas, font scale and Tol palette, so the figures sit together:
# nail magenta, mmseqs cyan and hmmer orange, as pid colours them. both hmmer
# arms are hmmer's orange, told apart by marker
SCALE = 1.75
mpl.rcParams.update({"font.size": mpl.rcParams["font.size"] * SCALE})
# 4:3, at pid's width
FIGSIZE = (12, 9)

ARMS = {
    "nail": "#EE3377",
    "mmseqs": "#33BBEE",
    "hmmer": "#EE7733",
    "hmmer-split": "#EE7733",
}
LABELS = {
    "nail": "nail",
    "mmseqs": "MMseqs2",
    "hmmer": "HMMER (one hmmsearch)",
    "hmmer-split": "HMMER (query split)",
}
MARKERS = {"nail": "o", "mmseqs": "o", "hmmer": "o", "hmmer-split": "v"}


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

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    for unit in units:
        threads_all = sorted({int(r["threads"]) for r in rows if r["unit"] == unit})
        series = []
        for arm in arms:
            pts = sorted(
                (int(r["threads"]), num(r, "wall_s"), num(r, "ideal_s"), num(r, "pct_ideal"))
                for r in rows
                if r["unit"] == unit and r["arm"] == arm and num(r, "wall_s")
            )
            if pts:
                series.append((arm, pts))

        suffix = f"-{unit}" if len(units) > 1 else ""
        wall(series, threads_all, out / f"thread-scaling-wall{suffix}.pdf")
        efficiency(series, threads_all, out / f"thread-scaling-efficiency{suffix}.pdf")


def style(arm):
    return dict(
        color=ARMS.get(arm, "#BBBBBB"),
        marker=MARKERS.get(arm, "o"),
        markersize=7,
        linewidth=2,
        label=LABELS.get(arm, arm),
    )


def frame(ax, threads_all):
    ax.set_xscale("log", base=2)
    ax.set_xticks(threads_all, [str(t) for t in threads_all])
    ax.minorticks_off()
    ax.set_xlabel("Threads")
    ax.grid(True)
    ax.set_axisbelow(True)


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print(f"wrote {path}")


def wall(series, threads_all, path):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    for arm, pts in series:
        t = [p[0] for p in pts]
        ax.plot(t, [p[1] for p in pts], **style(arm))
        # perfect scaling from this arm's lowest rung, in its colour
        ax.plot(t, [p[2] for p in pts], color=ARMS.get(arm, "#BBBBBB"), linestyle="--", linewidth=1)
    ax.plot([], [], color="#777777", linestyle="--", linewidth=1, label="perfect scaling")
    frame(ax, threads_all)
    ax.set_yscale("log")
    ax.set_ylabel("Wall time (s, log scale)")
    ax.legend(frameon=False)
    ax.set_title("Thread scaling: wall time")
    save(fig, path)


def efficiency(series, threads_all, path):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    for arm, pts in series:
        ax.plot([p[0] for p in pts], [p[3] for p in pts], **style(arm))
    ax.axhline(100.0, color="#999999", linestyle="--", linewidth=1)
    frame(ax, threads_all)
    ax.set_ylim(0, 105)
    ax.set_ylabel("Parallel efficiency (% of perfect speedup)")
    ax.legend(frameon=False, loc="lower left")
    ax.set_title("Thread scaling: parallel efficiency")
    save(fig, path)

if __name__ == "__main__":
    main()
