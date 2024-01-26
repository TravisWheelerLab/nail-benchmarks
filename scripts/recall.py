#! /bin/python3

import sys
import math
import subprocess
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# colors = [
#     "#1f78b4",  # (Blue)
#     "#33a02c",  # (Green)
#     "#e31a1c",  # (Red)
#     "#ff7f00",  # (Orange)
#     "#6a3d9a",  # (Purple)
#     "#a6cee3",  # (Light Blue)
#     "#b2df8a",  # (Light Green)
#     "#fb9a99",  # (Light Red)
#     "#fdbf6f",  # (Light Orange)
#     "#cab2d6",  # (Light Purple)
# ]

colors = [
    "#D81B60",  # red
    "#1E88E5",  # blue
    "#FFC107",  # yellow
    "#004D40",  # green
]

figsize = (10, 7)


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

        long_seq_target_paths = (
            benchmark_dir / "long-seq/target/").glob("*.fa")
        long_seq_query_paths = (
            benchmark_dir / "long-seq/query/").glob("*.fa")

        for (q, t) in zip(long_seq_query_paths, long_seq_target_paths):
            command = ["esl-seqstat", "-a", q]
            result = subprocess.run(
                command, stdout=subprocess.PIPE, text=True, check=True)
            lines = result.stdout.splitlines()
            for line in lines:
                if not line.startswith("="):
                    continue

                tokens = line.split()
                name = tokens[1]
                length = int(tokens[2])
                self.query_lengths[name] = length

            command = ["esl-seqstat", "-a", t]
            result = subprocess.run(
                command, stdout=subprocess.PIPE, text=True, check=True)
            lines = result.stdout.splitlines()
            for line in lines:
                if not line.startswith("="):
                    continue

                tokens = line.split()
                name = tokens[1]
                length = int(tokens[2])
                self.target_lengths[name] = length


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
        target_start=None,
        target_end=None,
        query_start=None,
        query_end=None,
        score=None,
        bias=None,
        cell_frac=None,
    ):
        self.target = target
        self.query = query
        self.target_start = target_start
        self.target_end = target_end
        self.query_start = query_start
        self.query_end = query_end
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
        target_start=None,
        target_end=None,
        query_start=None,
        query_end=None,
        score=None,
        bias=None,
        cell_frac=None,
    ):
        self.query_name = query_name
        self.target_name = target_name
        self.target_start = target_start,
        self.target_end = target_end,
        self.query_start = query_start,
        self.query_end = query_end,
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

                target_start = None
                if cols.target_start is not None:
                    target_start = int(line_tokens[cols.target_start])

                target_end = None
                if cols.target_end is not None:
                    target_end = int(line_tokens[cols.target_end])

                query_start = None
                if cols.query_start is not None:
                    query_start = int(line_tokens[cols.query_start])

                query_end = None
                if cols.query_end is not None:
                    query_end = int(line_tokens[cols.query_end])

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
                    target_start,
                    target_end,
                    query_start,
                    query_end,
                    score=score,
                    bias=bias,
                    cell_frac=cell_frac
                )

                if target.startswith("decoy"):
                    self.false_positives.append(hit)
                else:
                    target_tokens = target.split("/")
                    target_name = target_tokens[0]

                    if (target_start is not None and
                            target_end is not None and
                            len(target_tokens) == 3):
                        target_range = [int(i)
                                        for i in target_tokens[2].split('-')]
                        plant_start = target_range[0]
                        plant_end = target_range[1]
                        plant_length = plant_end - plant_start + 1

                        overlap_start = max(plant_start, target_start)
                        overlap_end = min(plant_end, target_end)

                        if overlap_start <= overlap_end:
                            overlap = overlap_end - overlap_start + 1
                        else:
                            overlap = 0

                        overlap_percentage = overlap / plant_length

                        # print(plant_start, plant_end)
                        # print(target_start, target_end)
                        # print(f"{overlap_percentage:3.2f}, {hit.evalue}")
                        # print()

                        if (target_name == query_name and
                                overlap_percentage >= 0.50):
                            self.true_positives.append(hit)
                        elif (target_name == query_name and
                                overlap_percentage == 0.0):

                            self.false_positives.append(hit)
                        else:
                            self.other.append(hit)

                    else:
                        if (target_name == query_name):
                            self.true_positives.append(hit)
                        else:
                            self.other.append(hit)

        def sort_key(h): return h.evalue

        self.true_positives.sort(key=sort_key)

        if (self.name == "hmmer" and
                cols.target_start is not None and
                cols.target_end is not None):
            unique_hits = {}
            for hit in self.true_positives:
                target_tokens = hit.target_name.split("/")
                target_name = target_tokens[0]

                if hit.query_name == target_name:
                    if ((hit.query_name, hit.target_name) not in unique_hits or
                            hit.evalue < unique_hits[(hit.query_name, hit.target_name)].evalue):
                        unique_hits[(hit.query_name, hit.target_name)] = hit

            self.true_positives = list(unique_hits.values())
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
    # cols = Cols(0, 2, 4)

    # best domain E-value
    # cols = Cols(0, 2, 7)
    paths = results_dir.glob("*.tbl")

    cols = Cols(0, 3, 12,
                target_start=17,
                target_end=18,
                query_start=15,
                query_end=16,
                score=7,
                )

    paths = results_dir.glob("*.domtbl")

    return [Hits(p, cols) for p in paths]


def read_mmseqs_results(results_dir):
    results_dir = results_dir / "mmseqs/"
    cols = Cols(0, 1, 6,
                target_start=2,
                target_end=3,
                query_start=4,
                query_end=5,
                )

    paths = results_dir.glob("*.tsv")

    return [Hits(p, cols) for p in paths]


def read_nail_results(results_dir):
    results_dir = results_dir / "nail/"
    cols = Cols(0, 1, 8,
                target_start=2,
                target_end=3,
                query_start=4,
                query_end=5,
                score=6,
                bias=7,
                cell_frac=9
                )

    paths = results_dir.glob("*.tsv")

    return [Hits(p, cols) for p in paths]


def plot_recall(hits, num_true_positives, num_queries):
    plt.close('all')
    plt.figure(figsize=figsize)

    hmmer_hits = next(filter(lambda h: h.name == "hmmer", hits))
    nail_hits = next(filter(lambda h: h.name == "nail default", hits))
    mmseqs_hits = next(filter(lambda h: h.name == "mmseqs default", hits))
    mmseqs_sensitive_hits = next(
        filter(lambda h: h.name == "mmseqs sensitive", hits))

    hits = [
        hmmer_hits,
        nail_hits,
        mmseqs_sensitive_hits,
        mmseqs_hits,
    ]

    labels = [
        "hmmer (default)",
        "nail (default)",
        "mmseqs (sensitive)",
        "mmseqs (default)",
    ]

    color = [
        colors[1],
        colors[2],
        colors[3],
        colors[0],
    ]

    ymin = 1.0
    ymax = 0.0
    for (h, c, l) in zip(hits, color, labels):
        (x, y, y_first, (x_fdr, y_fdr)) = h.recall_vs_mean_false(
            num_true_positives, num_queries)

        ymin = min(ymin, y_first)
        ymax = max(ymax, y[-1])

        plt.plot(
            x,
            y,
            color=c,
            label=l
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
            color=c,
        )

    plt.plot([], [], markeredgecolor='black', color='white', linestyle='', marker='D',
             label='1% False Discovery Rate')

    plt.plot([], [], color='black', linestyle='--',
             label='Recall Before First False Positive')

    ymin = (math.floor(ymin * 100) / 100) - 0.025
    ymax = (math.ceil(ymax * 100) / 100) + 0.025

    plt.xscale('log')

    plt.xlabel('Mean False Positives Per Search')
    plt.ylabel('Recall')
    plt.title('Pfam Domain Benchmark: Recall vs. Mean False Positives Per Search')

    plt.tick_params(
        axis='y',
        which='both',
        left=True,
        labelleft=True,
        right=True,
        labelright=True,
    )

    plt.xlim(1e-3, 1e1)
    plt.ylim(ymin, ymax)

    plt.legend(loc='upper left')

    plt.savefig("roc.pdf")
    # plt.show()


def plot_nail_bitscore(nail_hits):
    plt.close('all')
    plt.figure(figsize=figsize)

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

    n_greater = 0
    n_lower = 0
    losses = []
    gains = []
    both = []
    for (d, f) in zip(default_matched_hits, full_matched_hits):
        assert (d.target_name == f.target_name)

        perc = int((abs(f.score - d.score) / f.score) * 10000) / 100
        # perc = int((abs((f.score + f.bias) - (d.score + d.bias)) /
        #            (f.score + f.bias)) * 10000) / 100

        if d.score < f.score:
            n_lower += 1
            losses.append(perc)
        elif d.score > f.score:
            n_greater += 1
            gains.append(perc)

        both.append(perc)

    n = len(default_matched_hits)
    losses.sort()
    gains.sort()

    n_less_than_1_percent_loss = len(list(filter(lambda p: p < 1, losses)))
    n_less_than_2_percent_loss = len(list(filter(lambda p: p < 2, losses)))
    n_less_than_3_percent_loss = len(list(filter(lambda p: p < 3, losses)))
    n_less_than_5_percent_loss = len(list(filter(lambda p: p < 5, losses)))
    n_less_than_10_percent_loss = len(list(filter(lambda p: p < 10, losses)))

    n_less_than_1_percent_gain = len(list(filter(lambda p: p < 1, gains)))
    n_less_than_2_percent_gain = len(list(filter(lambda p: p < 2, gains)))
    n_less_than_3_percent_gain = len(list(filter(lambda p: p < 3, gains)))
    n_less_than_5_percent_gain = len(list(filter(lambda p: p < 5, gains)))
    n_less_than_10_percent_gain = len(list(filter(lambda p: p < 10, gains)))

    mean_loss = sum(losses) / len(losses)
    median_loss = losses[int(len(losses) / 2)]

    mean_gain = sum(gains) / len(gains)
    median_gain = gains[int(len(gains) / 2)]

    mean_score_sparse = sum([h.score for h in default_matched_hits]) / n
    mean_score_full = sum([h.score for h in full_matched_hits]) / n

    print(f" mean score (sparse): {mean_score_sparse:5.2f}")
    print(f" mean score (full): {mean_score_full:5.2f}")
    print()

    print(f"       n higher score: {n_greater} / {n}")
    print(f"  n less than 1% gain: {n_less_than_1_percent_gain} / {n_greater}")
    print(f"  n less than 2% gain: {n_less_than_2_percent_gain} / {n_greater}")
    print(f"  n less than 3% gain: {n_less_than_3_percent_gain} / {n_greater}")
    print(f"  n less than 5% gain: {n_less_than_5_percent_gain} / {n_greater}")
    print(
        f" n less than 10% gain: {n_less_than_10_percent_gain} / {n_greater}")
    print(f"      median increase: {median_gain:3.3f}%")
    print(f"        mean increase: {mean_gain:3.3f}%")
    print()
    print(f"        n lower score: {n_lower} / {n}")
    print(f"  n less than 1% loss: {n_less_than_1_percent_loss} / {n_lower}")
    print(f"  n less than 2% loss: {n_less_than_2_percent_loss} / {n_lower}")
    print(f"  n less than 3% loss: {n_less_than_3_percent_loss} / {n_lower}")
    print(f"  n less than 5% loss: {n_less_than_5_percent_loss} / {n_lower}")
    print(f" n less than 10% loss: {n_less_than_10_percent_loss} / {n_lower}")
    print(f"      median decrease: {median_loss:3.3f}%")
    print(f"        mean decrease: {mean_loss:3.3f}%")
    print()
    print(
        f"  n less than 1% diff: {n_less_than_1_percent_gain + n_less_than_1_percent_loss}")
    print(
        f"  n less than 2% diff: {n_less_than_2_percent_gain + n_less_than_2_percent_loss}")
    print(
        f"  n less than 3% diff: {n_less_than_3_percent_gain + n_less_than_3_percent_loss}")
    print(
        f"  n less than 5% diff: {n_less_than_5_percent_gain + n_less_than_5_percent_loss}")
    print(
        f" n less than 10% diff: {n_less_than_10_percent_gain + n_less_than_10_percent_loss}")

    # x = [h.score + h.bias for h in full_matched_hits]
    # y = [h.score + h.bias for h in default_matched_hits]

    x = [h.score for h in full_matched_hits]
    y = [h.score for h in default_matched_hits]

    max_x = max(x)
    max_y = max(y)
    max_val = max(max_x, max_y)

    coefficients = np.polyfit(x, y, deg=1)
    fit_line = np.poly1d(coefficients)

    plt.plot(x, fit_line(x), color=colors[2], label="Trend")
    plt.plot([0, max_val], [0, max_val], color=colors[0], label="y = x")

    plt.scatter(
        x,
        y,
        color=colors[1],
        marker='^',
        label='True Positives',
        s=10,
        alpha=0.8
    )

    plt.xlabel('Sequence Bitscore of Full Forward-Backward')
    plt.ylabel('Sequence Bitscore of Sparse Forward-Backward')

    # x = [h.score + h.bias for h in full_matched_hits]
    # y = [h.score + h.bias for h in default_matched_hits]

    x = [h.score for h in full_matched_hits]
    y = [h.score for h in default_matched_hits]

    max_x = max(x)
    max_y = max(y)
    max_val = max(max_x, max_y)

    coefficients = np.polyfit(x, y, deg=1)
    fit_line = np.poly1d(coefficients)

    plt.plot(x, fit_line(x), color=colors[2], label="Trend")
    plt.plot([0, max_val], [0, max_val], color=colors[0], label="y = x")

    plt.scatter(
        x,
        y,
        color=colors[1],
        marker='^',
        label='True Positives',
        s=10,
        alpha=0.8
    )

    plt.xlabel('Sequence Bitscore of Full Forward-Backward')
    plt.ylabel('Sequence Bitscore of Sparse Forward-Backward')
    plt.title('Pfam Domain Benchmark: Sequence Bitscore')
    plt.legend(loc='upper left')
    plt.grid()

    plt.tick_params(
        axis='y',
        which='both',
        left=True,
        labelleft=True,
        right=True,
        labelright=True,
    )

    plt.xlim([0, max_val])
    plt.ylim([0, max_val])

    plt.savefig("bitscore.pdf")
    # plt.show()


def plot_nail_cells(nail_hits, benchmark):
    plt.close('all')
    plt.figure(figsize=figsize)

    default_hits = next(
        filter(lambda h: h.name == "nail no-filters", nail_hits))
    long_seq_hits = next(filter(lambda h: h.name == "long-seq", nail_hits))

    hits_groups = [
        default_hits.false_positives,
        default_hits.true_positives,
        long_seq_hits.other,
    ]

    labels = [
        "Decoys",
        "True Positives",
        "Long Sequence",
    ]

    color = [
        colors[0],
        colors[1],
        colors[3],
    ]

    markers = [
        'v',
        '^',
        "*",
    ]

    alphas = [
        0.6,
        0.6,
        1.0,
    ]

    for (i, (hits, l, c, m, a)) in enumerate(zip(hits_groups, labels, color, markers, alphas)):
        x = []
        y = []
        for hit in hits:
            query_length = benchmark.query_lengths[hit.query_name]
            target_length = benchmark.target_lengths[hit.target_name]
            num_cells = query_length * target_length
            x.append(num_cells)
            y.append(hit.cell_frac)

        plt.scatter(
            x,
            y,
            color=c,
            marker=m,
            label=l,
            alpha=a
        )

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Total Cells Computed by Full Forward-Backward')
    plt.ylabel('Fraction of Cells Computed by Sparse Forward-Backward')
    plt.title('Pfam Domain Benchmark: Cells Computed')

    plt.tick_params(
        axis='y',
        which='both',
        left=True,
        labelleft=True,
        right=True,
        labelright=True,
    )

    plt.xlim(1e2, 1e10)
    plt.ylim(1e-5, 1.1)

    plt.legend(loc='upper right')
    plt.grid()

    plt.savefig("cells.png")


class Time:
    def __init__(self, path):
        self.name = path.name
        with open(path) as file:
            line = file.readline()
            tokens = line.split()
            self.seconds = float(tokens[1])


def plot_time(results_dir, hits, num_true_positives, num_queries):
    plt.close('all')
    plt.figure(figsize=figsize)

    hmmer_dir = results_dir / "hmmer/"
    hmmer_paths = hmmer_dir.glob("*.time")
    hmmer_times = {t.name: t for t in [Time(p) for p in hmmer_paths]}

    mmseqs_dir = results_dir / "mmseqs/"
    mmseqs_paths = mmseqs_dir.glob("*.time")
    mmseqs_times = {t.name: t for t in [Time(p) for p in mmseqs_paths]}

    nail_dir = results_dir / "nail/"
    nail_paths = nail_dir.glob("*.time")
    nail_times = {t.name: t for t in [Time(p) for p in nail_paths]}

    hmmer_time = hmmer_times["hmmer.8.time"].seconds,

    nail_default_time = nail_times["nail.seed.8.time"].seconds + \
        nail_times["nail.align.8.default.time"].seconds

    # nail_a8b12_time = nail_times["nail.seed.8.time"].seconds + \
    #     nail_times["nail.align.a8b12.time"].seconds

    nail_full_time = nail_times["nail.seed.8.time"].seconds + \
        nail_times["nail.align.full.time"].seconds

    mmseqs_default_time = mmseqs_times["mmseqs.default.8.time"].seconds

    mmseqs_sensitive_time = mmseqs_times["mmseqs.sensitive.8.time"].seconds

    mmseqs_nail_time = mmseqs_times["mmseqs.prefilter.nail.8.time"].seconds + \
        mmseqs_times["mmseqs.align.nail.8.time"].seconds

    times = [
        hmmer_time,
        nail_full_time,
        nail_default_time,
        # nail_a8b12_time,
        mmseqs_nail_time,
        mmseqs_sensitive_time,
        mmseqs_default_time,
    ]

    hmmer_hits = next(filter(lambda h: h.name == "hmmer", hits))
    nail_full_hits = next(filter(lambda h: h.name == "nail full", hits))
    nail_default_hits = next(filter(lambda h: h.name == "nail default", hits))
    # nail_a8b12_hits = next(filter(lambda h: h.name == "nail a8b12", hits))
    mmseqs_nail_hits = next(filter(lambda h: h.name == "mmseqs nail", hits))
    mmseqs_sensitive_hits = next(
        filter(lambda h: h.name == "mmseqs sensitive", hits))
    mmseqs_default_hits = next(
        filter(lambda h: h.name == "mmseqs default", hits))

    hits = [
        hmmer_hits,
        nail_full_hits,
        nail_default_hits,
        # nail_a8b12_hits,
        mmseqs_nail_hits,
        mmseqs_sensitive_hits,
        mmseqs_default_hits,
    ]

    recalls = [
        h.recall_vs_mean_false(
            num_true_positives, num_queries
        )[2] for h in hits
    ]

    labels = [
        "hmmsearch (default)",
        "nail (full DP)",
        "nail (default)",
        # "nail (alpha=12, beta=8)",
        "mmseqs (nail pipeline settings)",
        "mmseqs (sensitive)",
        "mmseqs (default)",
    ]

    color = [
        colors[1],
        colors[2],
        colors[2],
        colors[0],
        colors[0],
        colors[0],
    ]

    markers = [
        'o',
        'o',
        'D',
        'o',
        'D',
        's',
    ]

    for x, y, l, c, m in zip(recalls, times, labels, color, markers):
        plt.scatter(
            x,
            y,
            color=c,
            marker=m,
            label=l,
        )

    plt.xlabel('Recall before First False Positive')
    plt.ylabel('Runtime (sec)')
    plt.title('Pfam Domain Benchmark: Runtime vs Recall before First False Positive')

    plt.yscale('log')

    plt.tick_params(
        axis='y',
        which='both',
        left=True,
        labelleft=True,
        right=True,
        labelright=True,
    )

    plt.xlim([0.2, 0.8])
    plt.ylim([1e1, 10e3])

    plt.legend(loc='upper left')
    plt.grid()

    plt.savefig("runtime.pdf")
    # plt.show()


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

    plot_recall(
        all_hits, benchmark.num_true_positives, benchmark.num_queries)

    plot_time(results_dir, all_hits,
              benchmark.num_true_positives, benchmark.num_queries)

    plot_nail_bitscore(nail_hits)
    plot_nail_cells(nail_hits, benchmark)
