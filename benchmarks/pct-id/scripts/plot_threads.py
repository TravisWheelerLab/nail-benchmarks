#!/usr/bin/env python3

from plot import TOOL_COLORS

from collections import defaultdict
from pathlib import Path
import re
import os

import matplotlib.pyplot as plt


def scale_frac(threads, times, base):
    return [
        t / (base / th)
        for (th, t) in zip(threads, times)
    ]


def theoretical_time(threads, base):
    return [
        (base / t)
        for t in threads
    ]


def main(dir):
    curves = defaultdict(list)

    file_pat = re.compile(r"^(?P<params>.*)(\.prf|\.seq)\.time$")
    thread_pat = re.compile(r"t(?P<threads>\d+)")
    time_pat = re.compile(r"Elapsed.*\)\: (?P<time>.*)$")

    for path in os.listdir(dir):
        if m1 := file_pat.search(path):
            params = m1.group("params")
            tokens = params.split('-')

            threads = None
            for (i, tok) in enumerate(tokens):
                if m2 := thread_pat.search(tok):
                    threads = int(m2.group("threads"))
                    break

            tokens.pop(i)

            tool = "-".join(tokens)

            time = None
            with open(dir / path) as f:
                lines = f.readlines()
                times = []
                for line in lines:
                    if m := time_pat.search(line):
                        tokens = m.group("time").split(':')
                        tokens.reverse()
                        t = 0
                        for (i, tok) in enumerate(tokens):
                            t += float(tok) * (60.0 ** i)

                        times.append(t)
                time = max(times)

            curves[tool].append((threads, time))

    for tool in curves:
        curves[tool].sort(key=lambda x: x[0])

    curves = [(tool, *zip(*curves[tool])) for tool in curves]

#     for c in curves:
#         print(c)

    fig, ax = plt.subplots(figsize=(16, 9))
    plt.title("Runtime as threads increase")

    ax.set_xlabel("Threads")
    ax.set_xscale("log")

    ax.set_ylabel("Runtime (seconds)")
    ax.set_yscale("log")

    for label, x, y, in curves:
        # yy = scale_frac(x, y, y[0])
        yy = theoretical_time(x, y[0])

        print(label)

        for a, b in zip(y, yy):
            print(a, b)

        tool = label.split('-')[0]
        color = TOOL_COLORS[tool]
        ax.plot(
            y, yy,
            # color=color,
            marker="o",
            markersize=3,
        )

        # plt.text(
        #     x[-1] + 1, y[-1],
        #     label,
        #     color=color,
        #     fontweight="bold",
        #     ha="left",
        #     va="center",
        # )

    # fig.legend()

    plt.show()
    # plt.savefig(filename.with_suffix(".pdf"))


if __name__ == "__main__":
    import sys
    main(Path(sys.argv[1]))
