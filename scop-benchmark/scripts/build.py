import os
import argparse
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


decoy_cnt = 0


@dataclass
class Sequence:
    name: str
    seq: str

    def write(self, file, s, h, decoy=False):
        global decoy_cnt

        domain = s.domain_by_seq[self.name]
        family = h.family(domain)
        superfamily = h.superfamily(domain)
        fold = h.fold(domain)
        if decoy:
            file.write(f">decoy-{decoy_cnt}:{domain}|{family}|{superfamily}|{fold}\n{self.seq[::-1]}\n")
            decoy_cnt += 1
        else:
            file.write(f">{self.name}:{domain}|{family}|{superfamily}|{fold}\n{self.seq}\n")


class Domain:
    def __init__(self, tokens: [str]):
        self.sclass = tokens[3]
        self.fold = tokens[4]
        self.superfam = tokens[5]
        self.fam = tokens[6]
        self.name = tokens[0]

    def __eq__(self, o):
        return self.sclass == o.sclass and self.fold == o.fold and self.superfam == o.superfam and self.fam == o.fam and self.name == o.name

    def __str__(self):
        return f"{self.name} {self.fam} {self.superfam} {self.fold} {self.sclass}"


class Hierarchy:
    def __init__(self):
        self.domains = {}

    def add(self, tokens: [str]):
        domain = Domain(tokens)
        if domain.name in self.domains:
            assert (domain == self.domains[domain.name])
        else:
            self.domains[domain.name] = domain

    def family(self, domain: str):
        return self.domains[domain].fam

    def superfamily(self, domain: str):
        return self.domains[domain].superfam

    def fold(self, domain: str):
        return self.domains[domain].fold


class SequenceSet:
    def __init__(self, swipe_path, fasta_path, hierarchy):
        self.lengths_by_seq = {}
        self.domain_by_seq = {}
        self.seqs_by_domain = {}

        process = subprocess.Popen(
            ['esl-seqstat', '-a', fasta_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        for line in process.stdout:
            if line.startswith("="):
                tokens = line.split()
                self.lengths_by_seq[tokens[1]] = int(tokens[2])

        process.wait()

        domains_by_seq = {}
        with open(swipe_path, 'r') as file:
            for line in file:
                tokens = line.split()
                hierarchy.add(tokens)

                domain = tokens[0]
                seq = tokens[1]

                if seq in domains_by_seq:
                    domains_by_seq[seq].append(domain)
                else:
                    domains_by_seq[seq] = [domain]

        # filter for sequences that are only annotated by one domain
        for seq in domains_by_seq:
            domains = domains_by_seq[seq]

            if len(domains) == 1:
                domain = domains[0]
                if seq in self.domain_by_seq:
                    print("tried to set a domain to a seq twice")
                    exit(1)

                self.domain_by_seq[seq] = domain

                if domain in self.seqs_by_domain:
                    self.seqs_by_domain[domain].append(seq)
                else:
                    self.seqs_by_domain[domain] = [seq]

        self.seqs = [seq for seq in self.domain_by_seq]

    def domain_intersection(self, other):
        domains_a = set(self.seqs_by_domain)
        domains_b = set(other.seqs_by_domain)
        domains = list(domains_a & domains_b)
        return domains

    def sample(self, n):
        assert (n <= len(self.seqs))
        sample = np.random.choice(self.seqs, n, replace=False)
        return list(sample)

    def num_seqs_for(self, domain):
        return len(self.seqs_by_domain[domain])

    def sample_domain(self, domain, n):
        seqs = self.seqs_by_domain[domain]
        assert (n <= len(seqs))
        sample = np.random.choice(seqs, n, replace=False)
        return list(sample)


def fetch(seq_names, fasta_path, names_path="names.tmp"):
    with open(names_path, 'w') as file:
        for seq in seq_names:
            file.write(seq + '\n')

    process = subprocess.Popen(
        ['esl-sfetch', '-f', fasta_path, names_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    name = None
    seq = None
    seqs = []

    for line in process.stdout:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if seq is not None:
                seqs.append(Sequence(name, seq))
            name = line.lstrip(">").split()[0]
            seq = ""
        else:
            seq += line

    seqs.append(Sequence(name, seq))

    process.wait()

    return seqs


def length_histogram(lengths, title, bins='auto'):
    plt.hist(lengths, bins=bins, edgecolor='black')

    plt.xlabel('Length')
    plt.ylabel('Frequency')
    plt.title(f'{title} lengths')

    plt.savefig(f"{title}-lengths.pdf")


def length_stats(lengths):
    mean = np.mean(lengths)
    median = np.median(lengths)
    max_value = np.max(lengths)
    min_value = np.min(lengths)
    std_dev = np.std(lengths)

    print(f"mean:    {mean}")
    print(f"median:  {median}")
    print(f"max:     {max_value}")
    print(f"min:     {min_value}")
    print(f"std dev: {std_dev}")


if __name__ == "__main__":
    script_path = Path(os.path.abspath(__file__))
    data_path = script_path.parent / "../data/diamond-data/"

    query_path = data_path / "query_shuffled.faa"
    query_scop_path = data_path / "query_scop_annotation.tsv"
    target_path = data_path / "uniref50_annot_shuffled.faa"
    target_scop_path = data_path / "uniref50_scop_annotation.tsv"

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-s',
        '--random_seed',
        type=int,
        default=420,
        help='Random seed value'
    )

    parser.add_argument(
        '--num_domains',
        type=int,
        default=1_000,
        help='The number of query SCOP domains to use in the benchmark'
    )

    parser.add_argument(
        '--queries_per_domain',
        type=int,
        default=1,
        help='The number of sequences sampled per chosen query SCOP domain'
    )

    parser.add_argument(
        '--targets_per_query',
        type=int,
        default=10,
        help='The number of true match sequences sampled per query'
    )

    parser.add_argument(
        '--num_decoys',
        type=int,
        default=1_000_000,
        help='The total number of decoy sequences sampled'
    )

    parser.add_argument(
        '--name',
        type=str,
        default="benchmark",
        help='The benchmark name'
    )

    args = parser.parse_args()

    random_seed = args.random_seed
    np.random.seed(random_seed)

    benchmark_name = args.name
    benchmark_query_path = f"./{benchmark_name}.query.fa"
    benchmark_target_path = f"./{benchmark_name}.target.fa"

    os.makedirs(os.path.dirname(benchmark_query_path), exist_ok=True)

    if os.path.exists(benchmark_query_path):
        raise FileExistsError(f"The file '{benchmark_query_path}' already exists.")

    if os.path.exists(benchmark_target_path):
        raise FileExistsError(f"The file '{benchmark_target_path}' already exists.")

    hierarchy = Hierarchy()

    query_set = SequenceSet(query_scop_path, query_path, hierarchy)
    target_set = SequenceSet(target_scop_path, target_path, hierarchy)

    # only consider domains that appear in both the query & target set
    domains = query_set.domain_intersection(target_set)

    # only consider domains from which we can sample the desired number of targets
    domains = list(
        filter(
            lambda d: query_set.num_seqs_for(d) >= args.queries_per_domain
            and target_set.num_seqs_for(d) >= args.targets_per_query,
            domains
        )
    )

    # sort domains for determinism with random seed sample
    domains.sort()
    domains = np.random.choice(domains, args.num_domains, replace=False)

    queries = []
    targets = []
    for domain in domains:
        queries += query_set.sample_domain(domain, args.queries_per_domain)
        targets += target_set.sample_domain(domain, args.targets_per_query)

    decoys = target_set.sample(args.num_decoys)

    query_lengths = [query_set.lengths_by_seq[q] for q in queries]
    target_lengths = [target_set.lengths_by_seq[t] for t in targets]
    decoy_lengths = [target_set.lengths_by_seq[t] for t in decoys]

    # length_histogram(query_lengths, "query")
    # length_histogram(target_lengths, "target-true")
    # length_histogram(decoy_lengths, "target-decoy")

    with open(benchmark_query_path, "w") as file:
        seqs = fetch(queries, query_path)
        for seq in seqs:
            seq.write(file, query_set, hierarchy)

    with open(benchmark_target_path, "w") as file:
        seqs = fetch(targets, target_path)
        for seq in seqs:
            seq.write(file, target_set, hierarchy)

        seqs = fetch(decoys, target_path)
        for seq in seqs:
            seq.write(file, target_set, hierarchy, decoy=True)
