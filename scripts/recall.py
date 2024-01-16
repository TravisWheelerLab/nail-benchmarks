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


min_mean_false_positive = 1e-3
max_mean_false_positive = 10.0
num_mean_false_positive_points = 20

mean_false_positive_values = list(np.logspace(
    np.log10(min_mean_false_positive),
    np.log10(max_mean_false_positive),
    num_mean_false_positive_points
))

mean_false_positive_values.insert(0, 0)

recall_step = 0.01


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
        cell_frac=None,
    ):
        self.target = target
        self.query = query
        self.evalue = evalue
        self.score = score
        self.cell_frac = cell_frac


class Hit:
    def __init__(
        self,
        query_name,
        target_name,
        evalue,
        score=None,
        cell_frac=None,
    ):
        self.query_name = query_name
        self.target_name = target_name
        self.evalue = evalue
        self.score = score
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

                cell_frac = None
                if cols.cell_frac is not None:
                    cell_frac = float(line_tokens[cols.cell_frac])

                hit = Hit(
                    query_name,
                    target,
                    evalue,
                    score=score,
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
        return len(self.true_positives) + len(self.false_positives) + len(self.other)

    def recall_vs_mean_y_even(self, num_true_positives, num_queries):
        print(self.name)

        t = [h.target_name for h in self.true_positives]
        tp = list(set([h.target_name for h in self.true_positives]))

        assert (len(t) == len(tp))

        first_false_evalue = self.false_positives[0].evalue

        max_recall = len(self.true_positives) / num_true_positives

        true_before_first_false = list(
            filter(lambda x: x.evalue <= first_false_evalue, self.true_positives))

        num_true_before_first_false = len(true_before_first_false)
        recall_before_first_false = num_true_before_first_false / num_true_positives

        recall_start = math.ceil(recall_before_first_false * 100)
        recall_end = math.ceil(max_recall * 100)

        recall_values = [i / 100 for i in range(recall_start, recall_end)]

        recall_values.insert(0, recall_before_first_false)
        recall_values.append(max_recall)

        mean_false_values = []
        for r in recall_values:
            true_idx = (math.floor(r * num_true_positives) - 1)

            evalue_threshold = self.true_positives[true_idx].evalue
            false_below_thresh = list(
                filter(lambda x: x.evalue <= evalue_threshold, self.false_positives))

            num_false_below_thresh = len(false_below_thresh)
            mean_false_values.append(num_false_below_thresh / num_queries)

        x = mean_false_values
        y = recall_values

        return (x, y)

    def recall_vs_mean_false_x_even(self, num_true_positives, num_queries):
        # false indices maps to the point in the list of false positives
        # at which each approximate target mean false positive is achieved
        false_indices = []
        for target_value in mean_false_positive_values:
            if len(false_indices) > 0:
                start = false_indices[-1]
            else:
                start = 0

            for (idx, hit) in enumerate(self.false_positives[start:]):
                value = (idx + start) / num_queries
                if value >= target_value:
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
    cols = Cols(0, 1, 7, score=6, cell_frac=8)

    paths = results_dir.glob("*.tsv")

    return [Hits(p, cols) for p in paths]


def resample_points(x_vals, y_vals):
    x_resampled = []
    y_resampled = []
    for x_target in mean_false_positive_values:
        for (i, (x1, y1)) in enumerate(zip(x_vals, y_vals)):
            if x1 >= x_target:
                x0 = x_vals[i - 1]
                y0 = y_vals[i - 1]

                if x0 == x1:
                    y_target = y1
                    break
                else:
                    slope = (y1 - y0) / (x1 - x0)
                    y_target = y0 + (x_target - x0) * slope

                x_resampled.append(x_target)
                y_resampled.append(y_target)
                break

    return (x_resampled, y_resampled)


def plot_recall(hits, num_true_positives, num_queries):
    fig, ax = plt.subplots()

    ymin = 1.0
    ymax = 0.0
    for (h, c) in zip(hits, colors):
        vals = h.recall_vs_mean_y_even(num_true_positives, num_queries)
        (x, y) = resample_points(*vals)

        ymin = min(ymin, vals[1][0])
        ymax = max(ymax, vals[1][-1])

        ax.semilogx(
            *vals,
            color=c,
        )

        # ax.semilogx(
        #     x[1:],
        #     y[1:],
        #     marker='o',
        #     color=c,
        #     label=h.name
        # )

        ax.axhline(
            y=y[0],
            xmax=0,
            # linestyle='--',
            marker='D',
            color=c,
            clip_on=False
        )

    ymin = (math.floor(ymin * 100) / 100) - 0.025
    ymax = (math.ceil(ymax * 100) / 100) + 0.025

    ax.set_xlim([0, max_mean_false_positive])
    ax.set_ylim([ymin, ymax])

    ax.legend()

    if figures_path is not None:
        plt.savefig(figures_path)
    else:
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

    percentage_diffs = []
    for (d, f) in zip(default_matched_hits, full_matched_hits):
        assert (d.target_name == f.target_name)
        score_diff = abs(f.score - d.score)
        percentage_diff = score_diff / f.score
        percentage_diffs.append(percentage_diff)

    average_percentage_diff = sum(percentage_diffs) / len(percentage_diffs)

    x = [h.score for h in full_matched_hits]
    y = [h.score for h in default_matched_hits]

    plt.scatter(
        x,
        y,
        color='green',
        marker='^',
        label='True Positives',
        s=10,
        alpha=0.4
    )

    # Generating x values for the line
    x_line = np.linspace(min(x), max(x), 100)
    plt.plot(x_line, x_line, color='black',)

    plt.xlabel('Sequence Bitscore of Full Forward-Backward')
    plt.ylabel('Sequence Bitscore of Sparse Forward-Backward')
    plt.title('Pfam Domain Benchmark, Sequence Bitscore')
    plt.legend()
    plt.grid()
    plt.show()


def plot_nail_cells(nail_hits, benchmark):
    default_hits = next(filter(lambda h: h.name == "nail default", nail_hits))

    x_true_positive = []
    y_true_positive = []

    for hit in default_hits.true_positives:
        query_length = benchmark.query_lengths[hit.query_name]
        target_length = benchmark.target_lengths[hit.target_name]
        num_cells = query_length * target_length
        x_true_positive.append(num_cells)
        y_true_positive.append(hit.cell_frac)

    x_false_positive = []
    y_false_positive = []

    for hit in default_hits.false_positives:
        query_length = benchmark.query_lengths[hit.query_name]
        target_length = benchmark.target_lengths[hit.target_name]
        num_cells = query_length * target_length

        x_false_positive.append(num_cells)
        y_false_positive.append(hit.cell_frac)

    # fit lines
    x_values = [1e2, 1e9]

    true_coefficients = np.polyfit(
        np.log(x_true_positive), np.log(y_true_positive), deg=1)

    true_fit_line = np.exp(
        true_coefficients[1]) * x_values**true_coefficients[0]

    false_coefficients = np.polyfit(
        np.log(x_false_positive), np.log(y_false_positive), deg=1)

    false_fit_line = np.exp(
        false_coefficients[1]) * x_values**false_coefficients[0]

    # plotting

    plt.loglog(x_values, false_fit_line, color='red')

    plt.loglog(x_values, true_fit_line, color='green')

    plt.loglog(
        x_false_positive,
        y_false_positive,
        color='red',
        marker='^',
        label='False Positives',
        linestyle='',
        markersize=5,
        alpha=0.4
    )

    plt.loglog(
        x_true_positive,
        y_true_positive,
        color='green',
        marker='^',
        label='True Positives',
        linestyle='',
        markersize=5,
        alpha=0.4
    )

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

    # plot_recall(all_hits, benchmark.num_true_positives, benchmark.num_queries)
    # plot_nail_bitscore(nail_hits)
    plot_nail_cells(nail_hits, benchmark)
