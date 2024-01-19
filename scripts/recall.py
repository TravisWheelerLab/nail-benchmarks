#! /bin/python3

import sys
import math
import subprocess
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

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


class Benchmark:
    def __init__(self, benchmark_dir):
        benchmark_name = benchmark_dir.name

        # positive targets
        benchmark_pos = benchmark_dir / f"{benchmark_name}.pos"

        self.positives = []
        with open(benchmark_pos) as file:
            for line in file:
                positive = Positive(line)
                self.positives.append(positive)

        self.queries = list(set([p.name for p in self.positives]))

        self.num_queries = len(self.queries)
        self.num_true_positives = len(self.positives)

        # target lengths
        benchmark_target_fa = benchmark_dir / f"{benchmark_name}.test.fa"
        command = ["esl-seqstat", "-a", benchmark_target_fa]
        result = subprocess.run(
            command, stdout=subprocess.PIPE, text=True, check=True)

        lines = result.stdout.splitlines()

        self.target_lengths = {}
        for line in lines:
            if not line.startswith("="):
                continue

            tokens = line.split()
            name = tokens[1]
            length = int(tokens[2])
            self.target_lengths[name] = length

        # query (model) lengths
        benchmark_query_hmm = benchmark_dir / f"{benchmark_name}.train.hmm"

        command = ["hmmstat", benchmark_query_hmm]
        result = subprocess.run(
            command, stdout=subprocess.PIPE, text=True, check=True)

        lines = result.stdout.splitlines()

        self.query_lengths = {}
        for line in lines:
            if line.startswith("#"):
                continue

            tokens = line.split()
            name = tokens[1]
            length = int(tokens[3])
            self.query_lengths[name] = length


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
        score=None,
        bias=None,
        cell_frac=None,
    ):
        self.target = target
        self.query = query
        self.evalue = evalue
        self.score = score
        self.bias = bias
        self.cell_frac = cell_frac


class Hit:
    def __init__(
        self,
        query_name,
        target_name,
        evalue,
        score=None,
        bias=None,
        cell_frac=None,
    ):
        self.query_name = query_name
        self.target_name = target_name
        self.evalue = evalue
        self.score = score
        self.bias = bias
        self.cell_frac = cell_frac

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

                score = None
                if cols.score is not None:
                    score = float(line_tokens[cols.score])

                bias = None
                if cols.bias is not None:
                    bias = float(line_tokens[cols.bias])

                cell_frac = None
                if cols.cell_frac is not None:
                    cell_frac = float(line_tokens[cols.cell_frac])

                hit = Hit(
                    query_name,
                    target,
                    evalue,
                    score=score,
                    bias=bias,
                    cell_frac=cell_frac
                )

                if target.startswith("decoy"):
                    self.false_positives.append(hit)
                else:
                    target_tokens = target.split("/")
                    target_name = target_tokens[0]

                    if (target_name == query_name):
                        self.true_positives.append(hit)
                    else:
                        self.other.append(hit)

        def sort_key(h): return h.evalue

        self.true_positives.sort(key=sort_key)
        self.false_positives.sort(key=sort_key)
        self.other.sort(key=sort_key)

    def num_hits(self):
        return (
            len(self.true_positives)
            + len(self.false_positives)
            + len(self.other)
        )

    def recall_vs_mean_false(self, num_true_positives, num_queries):
        t = [h.target_name for h in self.true_positives]
        tp = list(set([h.target_name for h in self.true_positives]))

        assert (len(t) == len(tp))

        x = []
        y = []

        combined_hits = [(h, True) for h in self.true_positives] + \
            [(h, False) for h in self.false_positives]

        combined_hits.sort(key=lambda h: h[0].evalue)

        true_count = 0
        false_count = 0
        fdr_point = None

        for (hit, is_true) in combined_hits:
            if is_true:
                true_count += 1
            else:
                false_count += 1

            recall = true_count / num_true_positives
            mean_false_positive = false_count / num_queries

            x.append(mean_false_positive)
            y.append(recall)

            if fdr_point is None:
                false_discovery_rate = false_count / true_count
                if abs(0.01 - false_discovery_rate) < 1e-3:
                    fdr_point = (mean_false_positive, recall)

        for (idx, val) in enumerate(x):
            if val > 0:
                y_first = y[idx - 1]
                break

        return (x, y, y_first, fdr_point)


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
    cols = Cols(0, 1, 8, score=6, bias=7, cell_frac=9)

    paths = results_dir.glob("*.tsv")

    return [Hits(p, cols) for p in paths]


def plot_recall(hits, num_true_positives, num_queries):
    plotted = [
        "hmmer",
        "nail default",
        "mmseqs sensitive",
        "mmseqs default",
    ]

    hits = list(filter(lambda h: h.name in plotted, hits))

    ymin = 1.0
    ymax = 0.0
    for (h, c) in zip(hits, colors):
        (x, y, y_first, (x_fdr, y_fdr)) = h.recall_vs_mean_false(
            num_true_positives, num_queries)

        ymin = min(ymin, y_first)
        ymax = max(ymax, y[-1])

        plt.plot(
            x,
            y,
            color=c,
            label=h.name
        )

        plt.plot(
            x_fdr,
            y_fdr,
            color=c,
            marker='D',
            markeredgecolor='black',
        )

        plt.axhline(
            y=y_first,
            xmax=10,
            linestyle='--',
            alpha=0.4,
            # marker='D',
            color=c,
            # clip_on=False
        )

    plt.plot([], [], color='black', linestyle='', marker='D',
             label='1% False Discovery Rate')

    plt.plot([], [], color='black', linestyle='--',
             label='Recall Before First False Positive')

    ymin = (math.floor(ymin * 100) / 100) - 0.025
    ymax = (math.ceil(ymax * 100) / 100) + 0.025

    plt.xscale('log')

    plt.xlabel('Mean False Positives Per Search')
    plt.ylabel('Recall')
    plt.title('Pfam Domain Benchmark, Recall vs. Mean False Positives Per Search')

    plt.xlim(1e-3, 1e1)
    plt.ylim(ymin, ymax)

    plt.legend()

    plt.savefig("recall-vs-mean-false.pdf")
    plt.show()


def plot_nail_bitscore(nail_hits):
    default_hits = next(filter(lambda h: h.name == "nail default", nail_hits))
    full_hits = next(filter(lambda h: h.name == "nail full", nail_hits))

    # take the target names of all of the default true positives
    target_names = [h.target_name for h in default_hits.true_positives]

    # find all of the full DP hits that have the same target name
    full_matched_hits = list(filter(
        lambda h: h.target_name in target_names,
        full_hits.true_positives)
    )

    # we should never have more full DP hits matched to the default hits
    assert (len(full_matched_hits) <= len(default_hits.true_positives))

    # just in case there were some default hits that don't match to a full
    # DP hit, take the set of target names from the matched full hits
    full_matched_target_names = [h.target_name for h in full_matched_hits]

    # now subset the default hits using the reduced target names
    default_matched_hits = list(filter(
        lambda h: h.target_name in full_matched_target_names,
        default_hits.true_positives)
    )

    # now we should have the same number of default & full DP hits
    assert (len(default_matched_hits) == len(full_matched_hits))

    def s(h): return h.target_name
    default_matched_hits.sort(key=s)
    full_matched_hits.sort(key=s)

    for (d, f) in zip(default_matched_hits, full_matched_hits):
        assert (d.target_name == f.target_name)

    # x = [h.score + h.bias for h in full_matched_hits]
    # y = [h.score + h.bias for h in default_matched_hits]

    x = [h.score for h in full_matched_hits]
    y = [h.score for h in default_matched_hits]

    max_x = max(x)
    max_y = max(y)
    max_val = max(max_x, max_y)

    coefficients = np.polyfit(x, y, deg=1)
    fit_line = np.poly1d(coefficients)

    plt.plot(x, fit_line(x), color=colors[3], label="Trend")
    plt.plot([0, max_val], [0, max_val], color=colors[1], label="y = x")

    plt.scatter(
        x,
        y,
        color=colors[0],
        marker='^',
        label='True Positives',
        s=10,
        alpha=0.4
    )

    plt.xlabel('Sequence Bitscore of Full Forward-Backward')
    plt.ylabel('Sequence Bitscore of Sparse Forward-Backward')
    plt.title('Pfam Domain Benchmark, Sequence Bitscore')
    plt.legend()
    plt.grid()
    plt.show()


def plot_nail_cells(nail_hits, benchmark):
    default_hits = next(filter(lambda h: h.name == "nail default", nail_hits))

    hits_groups = [
        default_hits.false_positives,
        default_hits.true_positives,
    ]

    labels = [
        "Decoys",
        "True Positives",
    ]

    color = [
        colors[2],
        colors[0],
    ]

    for (hits, l, c) in zip(hits_groups, labels, color):
        x = []
        y = []
        for hit in hits:
            query_length = benchmark.query_lengths[hit.query_name]
            target_length = benchmark.target_lengths[hit.target_name]
            num_cells = query_length * target_length
            x.append(num_cells)
            y.append(hit.cell_frac)

        coefficients = np.polyfit(np.log10(x), np.log10(y), deg=1)

        x_fit = [1e2, 1e9]
        y_fit = np.exp(coefficients[1]) * x_fit**coefficients[0]

        plt.plot(x_fit, y_fit, color=c, zorder=0)

        plt.scatter(
            x,
            y,
            color=c,
            marker='^',
            label=l,
            linestyle='',
            alpha=0.4
        )

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Total Cells Computed by Full Forward-Backward')
    plt.ylabel('Fraction of Cells Computed by Sparse Forward-Backward')
    plt.title('Pfam Domain Benchmark, Cells Computed')

    plt.xlim(1e2, 1e9)
    plt.ylim(1e-5, 1)

    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    figures_path = None
    if len(sys.argv) < 2:
        print("usage: ./recall.py <benchmark_dir> [figures/]")
        exit()
    elif len(sys.argv) > 2:
        figures_path = sys.argv[2]

    benchmark_dir = Path(sys.argv[1])
    results_dir = benchmark_dir / "results/"

    hmmer_hits = read_hmmer_results(results_dir)
    nail_hits = read_nail_results(results_dir)
    mmseqs_hits = read_mmseqs_results(results_dir)

    all_hits = hmmer_hits + nail_hits + mmseqs_hits

    benchmark = Benchmark(benchmark_dir)

    plot_recall(all_hits, benchmark.num_true_positives, benchmark.num_queries)
    plot_nail_bitscore(nail_hits)
    plot_nail_cells(nail_hits, benchmark)
