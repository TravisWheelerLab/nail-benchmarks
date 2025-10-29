#! /usr/bin/python3

from common import FastaIndex

import sys
import os


def write(fa_path, out_dir, splits):
    os.makedirs(out_dir, exist_ok=True)
    _, ext = os.path.splitext(fa_path)

    fa_lines = [line for line in open(fa_path)]

    for (spl_idx, spl) in enumerate(splits):
        with open(f"{out_dir}/{spl_idx}{ext}", "w") as f:
            for (start, end, _) in spl:
                for line in fa_lines[start:end + 1]:
                    f.write(line)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: fastabalance <file> <num_splits> [out_dir]")
        sys.exit(1)

    fa_path = sys.argv[1]
    n = int(sys.argv[2])
    index = FastaIndex(fa_path)

    index.ranges.sort(key=lambda x: x[2])
    splits = index.split(n)

    if len(sys.argv) > 3:
        out_dir = sys.argv[3]
    else:
        name, _ = os.path.splitext(fa_path)
        out_dir = f"{name}-splits"

    write(fa_path, out_dir, splits)
