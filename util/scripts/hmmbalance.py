#! /usr/bin/python3

from common import HmmIndex

import sys
import os


def write(hmm_path, out_dir, splits):
    os.makedirs(out_dir, exist_ok=True)
    _, ext = os.path.splitext(hmm_path)

    hmm_lines = [line for line in open(hmm_path)]

    for (spl_idx, spl) in enumerate(splits):
        with open(f"{out_dir}/{spl_idx}{ext}", "w") as f:
            for (start, end, _, _) in spl:
                for line in hmm_lines[start:end + 1]:
                    f.write(line)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: hmmbalance <file> <num_splits> [out_dir]")
        sys.exit(1)

    hmm_path = sys.argv[1]
    n = int(sys.argv[2])
    index = HmmIndex(hmm_path)

    index.ranges.sort(key=lambda x: x[2])
    splits = index.split(n)

    if len(sys.argv) > 3:
        out_dir = sys.argv[3]
    else:
        name, _ = os.path.splitext(hmm_path)
        out_dir = f"{name}-splits"

    write(hmm_path, out_dir, splits)
