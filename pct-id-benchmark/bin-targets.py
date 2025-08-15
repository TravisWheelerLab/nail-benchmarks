#! /usr/bin/python3

import sys
import os
import subprocess
from pathlib import Path
from collections import defaultdict

from boio import Fasta, Stockholm


AMINO_ALPH = "ACDEFGHIKLMNPQRSTVWYXB"


def compute_pid_no_indels(s1: str, s2: str) -> float:
    assert (len(s1) == len(s2))
    s1 = s1.upper()
    s2 = s2.upper()

    match_cnt = 0
    pos_cnt = 0

    for (a, b) in zip(s1, s2):
        if a in AMINO_ALPH and b in AMINO_ALPH:
            pos_cnt += 1
            match_cnt += a == b

    return float(match_cnt) / float(pos_cnt)


def compute_pid_with_indels(s1: str, s2: str) -> float:
    assert (len(s1) == len(s2))
    s1 = s1.upper()
    s2 = s2.upper()

    match_cnt = 0
    pos_cnt = 0

    for (a, b) in zip(s1, s2):
        if a in AMINO_ALPH or b in AMINO_ALPH:
            pos_cnt += 1
            if a == b:
                match_cnt += 1

    return float(match_cnt) / float(pos_cnt)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("usage: python bin_targets.py <test.fa> <train.sto> <source.sto> <out/>")
        sys.exit(1)

    test_fa_path = Path(sys.argv[1])
    query_sto_path = Path(sys.argv[2])
    src_sto_path = Path(sys.argv[3])
    out_dir = Path(sys.argv[4])

    test_fa = Fasta.from_path(test_fa_path)
    query_sto = Stockholm.parse(query_sto_path)
    src_sto = Stockholm.parse(src_sto_path)

    test_seqs_by_query = defaultdict(list)

    for t_seq in test_fa:
        q_fam_name = t_seq.name.split("/")[0]
        test_seqs_by_query[q_fam_name].append(t_seq)

    best_pid_by_target = {}

    for q_fam_name in query_sto.records:
        names_in_query = [s for s in query_sto.records[q_fam_name].sequences]
        names_in_target = [s.extra.split()[-1] for s in test_seqs_by_query[q_fam_name]]

        seqs = {
            k: v.sequence
            for k, v in src_sto.records[q_fam_name].sequences.items()
            if k in names_in_query or k in names_in_target
        }

        for t_name in names_in_target:
            t_seq = seqs[t_name]
            best = (0.0, "")
            for query_seq_name in names_in_query:
                q_seq = seqs[query_seq_name]
                pid = compute_pid_with_indels(q_seq, t_seq)
                if pid > best[0]:
                    best = (pid, query_seq_name)

            best_pid_by_target[t_name] = best

    fams_dir = out_dir / "fams/"
    os.makedirs(fams_dir, exist_ok=True)

    for (q_fam_name, seqs) in test_seqs_by_query.items():
        query_dir = fams_dir / f"{q_fam_name}"

        os.makedirs(query_dir, exist_ok=True)

        for i, t_seq in enumerate(seqs):
            q_sto = query_sto.records[q_fam_name]
            domain = t_seq.extra.split()[-1]
            (pid, query_seq_name) = best_pid_by_target[domain]
            q_seq = q_sto.sequences[query_seq_name].fasta_record()

            q_seq.extra += f" | pid: {pid:.2f} | to: {t_seq.name}"
            t_seq.extra += f" | pid: {pid:.2f} | to: {query_seq_name}"

            sto_path = query_dir / "query.sto"
            with open(sto_path, "w") as f:
                f.write(str(q_sto))

            hmm_path = query_dir / "query.hmm"
            subprocess.run(
                ["hmmbuild", hmm_path, sto_path],
                check=True,
                stdout=subprocess.DEVNULL,
            )

            cons_path = query_dir / "consensus.fa"
            with open(cons_path, "w") as f:
                subprocess.run(
                    ["hmmemit", "-c", hmm_path],
                    check=True,
                    stdout=f,
                )

            with open(query_dir / f"{i}.q.fa", "w") as f:
                f.write(str(q_seq))

            with open(query_dir / f"{i}.t.fa", "w") as f:
                f.write(str(t_seq))
