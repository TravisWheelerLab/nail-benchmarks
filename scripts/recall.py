#! /bin/python3

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FuncFormatter

colors = [
    "#1f78b4",  # (Blue)
    "#33a02c",  # (Green)
    "#e31a1c",  # (Red)
    "#ff7f00",  # (Orange)
    "#6a3d9a",  # (Purple)
    "#a6cee3",  # (Light Blue)
    "#b2df8a",  # (Light Green)
    "#fb9a99",  # (Light Red)
    "#fdbf6f",  # (Light Orange)
    "#cab2d6",  # (Light Purple)
]


min_mean_false_positive = 1e-3
max_mean_false_positive = 10.0
num_mean_false_positive_points = 20

mean_false_positive_values = np.logspace(
    np.log10(min_mean_false_positive),
    np.log10(max_mean_false_positive),
    num_mean_false_positive_points
)

recall_step = 0.01


class Positive:
    def __init__(self, line):
        tokens = line.split()
        self.label = tokens[0]

        label_tokens = tokens[0].split("/")
        self.name = label_tokens[0]
        self.id = int(label_tokens[1])

        coord_tokens = label_tokens[2].split("-")
        self.start = int(coord_tokens[0])
        self.end = int(coord_tokens[1])

    def __str__(self):
        return f"{self.name}\t{self.id}\t{self.start}..{self.end}"


class Cols:
    def __init__(
        self,
        target,
        query,
        evalue,
    ):
        self.target = target
        self.query = query
        self.evalue = evalue


class Hit:
    def __init__(
        self,
        query_name,
        target_name,
        evalue,
    ):
        self.query_name = query_name
        self.target_name = target_name
        self.evalue = evalue

    def __str__(self):
        return f"{self.query_name} {self.target_name} {self.evalue}"


class Hits:
    def __init__(self, path, cols):
        self.name = " ".join(path.name.split(".")[:-1])
        self.true_positives = []
        self.false_positives = []
        self.other = []

        with open(path) as file:
            for line in file:
                if line.startswith("#"):
                    continue

                line_tokens = line.split()
                # the target name is formatted like:
                #     <target_name>/<id>/<from>-<to>, or
                #     <decoy#>
                target = line_tokens[cols.target]
                query_name = line_tokens[cols.query]

                evalue = float(line_tokens[cols.evalue])

                if target.startswith("decoy"):
                    self.false_positives.append(
                        Hit(query_name, target, evalue)
                    )
                else:
                    target_tokens = target.split("/")
                    target_name = target_tokens[0]
                    # target_id = int(target_tokens[1])

                    hit = Hit(query_name, target, evalue)
                    if (target_name == query_name):
                        self.true_positives.append(hit)
                    else:
                        self.other.append(hit)

        def sort_key(h): return h.evalue

        self.true_positives.sort(key=sort_key)
        self.false_positives.sort(key=sort_key)
        self.other.sort(key=sort_key)

    def recall_vs_mean_false(self, num_true_positives, num_queries):
        # false indices maps to the point in the list of false positives
        # at which each approximate target mean false positive is achieved
        false_indices = [0]
        for target_value in mean_false_positive_values:
            start = false_indices[-1]

            for (idx, hit) in enumerate(self.false_positives[start:]):
                value = idx/num_queries
                if value >= target_value:
                    # print("{:.1e}\t{:.1e}\t{:.1e}".format(
                    #     target_value,
                    #     value,
                    #     abs(value - target_value)
                    # ))
                    false_indices.append(idx + start)
                    break

        evalue_thresholds = [
            self.false_positives[idx].evalue for idx in false_indices
        ]

        true_indices = [0]
        for e in evalue_thresholds:
            start = true_indices[-1]
            for (idx, hit) in enumerate(self.true_positives[start:]):
                if hit.evalue >= e:
                    true_indices.append(idx + start)
                    break

        true_indices = true_indices[1:]

        recall_values = [
            idx / num_true_positives for idx in true_indices][:20]

        y_len = len(recall_values)
        x = mean_false_positive_values[:y_len]
        y = recall_values
        return (x, y)

        pass

    def __str__(self):
        return f"TP:{len(self.true_positives)}"
        + f"FP:{len(self.false_positives)}"
        + f"O:{len(self.other)}"


def read_hmmer_results(results_dir):
    results_dir = results_dir / "hmmer/"
    # full seq E-value
    cols = Cols(0, 2, 4)

    # best domain E-value
    # cols = Cols(0, 2, 7)

    paths = results_dir.glob("*.tbl")

    return [Hits(p, cols) for p in paths]


def read_mmseqs_results(results_dir):
    results_dir = results_dir / "mmseqs/"
    cols = Cols(0, 1, 6)

    paths = results_dir.glob("*.tsv")

    return [Hits(p, cols) for p in paths]


def read_nail_results(results_dir):
    results_dir = results_dir / "nail/"
    cols = Cols(0, 1, 7)

    paths = results_dir.glob("*.tsv")

    return [Hits(p, cols) for p in paths]


def read_benchmark_info(benchmark_dir):
    benchmark_name = benchmark_dir.name
    benchmark_pos = benchmark_dir / f"{benchmark_name}.pos"

    positives = []
    with open(benchmark_pos) as file:
        for line in file:
            positive = Positive(line)
            positives.append(positive)

    queries = list(set([p.name for p in positives]))

    num_queries = len(queries)
    num_true_positives = len(positives)
    return (num_queries, num_true_positives)


def plot():
    if len(sys.argv) < 2:
        print("usage: ./recall.py <benchmark_dir>")
        exit()

    benchmark_dir = Path(sys.argv[1])
    results_dir = benchmark_dir / "results/"

    hits = []

    hits += read_hmmer_results(results_dir)
    hits += read_nail_results(results_dir)
    hits += read_mmseqs_results(results_dir)

    (num_queries, num_true_positives) = read_benchmark_info(benchmark_dir)

    for (h, c) in zip(hits, colors):
        x = h.recall_vs_mean_false(num_true_positives, num_queries)
        plt.semilogx(*x, 'o', color=c, label=h.name)
        plt.semilogx(*x, color=c)

    # Format x-axis ticks as floats
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    plt.gca().xaxis.set_major_formatter(
        FuncFormatter(lambda x, _: '{:.5g}'.format(x)))
    plt.gca().yaxis.set_major_formatter(
        FuncFormatter(lambda y, _: '{:.5g}'.format(y)))

    plt.xlim([0.0, max_mean_false_positive])

    ymin = 0.1
    ymax = 0.9
    plt.ylim([ymin, ymax])

    plt.legend()
    plt.show()
    # plt.savefig("recall-vs-mean-false-positive.png")


if __name__ == "__main__":
    plot()
