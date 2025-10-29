#! /usr/bin/python3

from common import FastaIndex

import argparse
import random
import time

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("fa_path")
    p.add_argument("-n", type=int, required=True)
    p.add_argument("--seed", type=int)
    args = p.parse_args()

    if args.seed is not None:
        seed = args.seed
    else:
        seed = int(time.time_ns() & 0xFFFFFFFF)

    random.seed(args.seed)

    index = FastaIndex(args.fa_path)
    ranges = [index.ranges[n] for n in random.sample(range(0, len(index.ranges)), args.n)]

    for r in ranges:
        print(index.read_lines(r[0], r[1]), end="")
