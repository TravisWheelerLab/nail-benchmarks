#! /usr/bin/python3

from common import HmmIndex

import argparse
import os
from pathlib import Path

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("hmm_path")
    p.add_argument("--min", type=int, default=50)
    p.add_argument("--max", type=int, default=1000)
    p.add_argument("--step", type=int, default=50)
    p.add_argument("--out", type=str)
    args = p.parse_args()

    index = HmmIndex(args.hmm_path)

    ranges_by_length = {}
    for r in index.ranges:
        l = r[2]
        if l in ranges_by_length:
            ranges_by_length[l].append(r)
        else:
            ranges_by_length[l] = [r]

    if args.out is not None:
        out_dir = Path(f"{args.out}")
    else:
        out_dir = Path(f"{args.hmm_path}-bins/")
    out = os.makedirs(out_dir, exist_ok=True)

    for upper in range(args.min, args.max + args.step, args.step):
        lower = upper - args.step + 1
        with open(out_dir / f"{lower}-{upper}.hmm", "w") as out:
            for l in range(lower, upper + 1):
                if l in ranges_by_length:
                    for r in ranges_by_length[l]:
                        out.write(index.read_lines(r[0], r[1]))
