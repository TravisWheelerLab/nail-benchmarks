#! /usr/bin/python3

from common import HmmIndex

import argparse
import math
import os
import random
from pathlib import Path


def bindex(length):
    x = length % args.step
    y = length - x
    return int(y / args.step) - int(x == 0)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("hmm_path")
    p.add_argument("--max", type=int, default=1000)
    p.add_argument("--step", type=int, default=50)
    p.add_argument("--take", type=int)
    p.add_argument("--out", type=str)
    args = p.parse_args()

    index = HmmIndex(args.hmm_path)

    bins = [(upper - args.step + 1, upper) for upper in range(args.step, args.max + args.step, args.step)]
    n_bins = len(bins)
    ranges_by_bin = [[] for _ in range(n_bins)]
    for r in index.ranges:
        length = r[2]

        if length > args.max:
            continue

        ranges_by_bin[bindex(length)].append(r)

    if args.out is not None:
        out_dir = Path(f"{args.out}")
    else:
        out_dir = Path(f"{args.hmm_path}-bins/")
    out = os.makedirs(out_dir, exist_ok=True)

    if args.take is not None:
        take = args.take
    else:
        take = math.inf

    for (b, bin) in enumerate(bins):
        sz = len(ranges_by_bin[b])
        if sz > take:
            sample = [ranges_by_bin[b][idx] for idx in random.sample(range(0, sz), take)]
        else:
            sample = ranges_by_bin[b]

        with open(out_dir / f"{bin[0]}-{bin[1]}.hmm", "w") as out:
            for r in sample:
                out.write(index.read_lines(r[0], r[1]))
