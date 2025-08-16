#! /usr/bin/python3

import sys
import os
import subprocess
from pathlib import Path
from collections import defaultdict

from boio import Fasta, Stockholm


AMINO_ALPH = "ACDEFGHIKLMNPQRSTVWYXB"


def compute_pid(s1: str, s2: str) -> float:
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
    if len(sys.argv) != 2:
        print("usage: python bin_targets.py <benchmark/>")
        sys.exit(1)

    dir = Path(sys.argv[1])
    pm_dir = dir / "profmark"
    target_fa_path = pm_dir / "target.fa"
    query_sto_path = pm_dir / "query.sto"
    src_sto_path = pm_dir / "source.sto"
    out_dir = dir / "fams/"

    os.makedirs(out_dir, exist_ok=True)

    target_fa = Fasta.from_path(target_fa_path)
    query_sto = Stockholm.from_path(query_sto_path)
    src_sto = Stockholm.from_path(src_sto_path)

    target_seqs_by_fam = defaultdict(list)

    for t_seq in target_fa:
        fam = t_seq.name.split("/")[0]
        target_seqs_by_fam[fam].append(t_seq)
    

    best_pid_by_target = {}

    for i, fam in enumerate(query_sto.records):
        names_in_query = [s for s in query_sto.records[fam].sequences]
        names_in_target = [s.extra.split()[-1] for s in target_seqs_by_fam[fam]]

        seqs = {
            k: v.sequence
            for k, v in src_sto.records[fam].sequences.items()
            if k in names_in_query or k in names_in_target
        }

        for t_name in names_in_target:
            t_seq = seqs[t_name]
            best = (0.0, "")
            for query_seq_name in names_in_query:
                q_seq = seqs[query_seq_name]
                pid = compute_pid(q_seq, t_seq)
                if pid > best[0]:
                    best = (pid, query_seq_name)

            best_pid_by_target[t_name] = best

    exit()

    for (fam, seqs) in target_seqs_by_fam.items():
        fam_dir = out_dir / f"{fam}"
        os.makedirs(fam_dir, exist_ok=True)

        q_sto = query_sto.records[fam]

        sto_path = fam_dir / "query.sto"
        with open(sto_path, "w") as f:
            f.write(str(q_sto))

        hmm_path = fam_dir / "query.hmm"
        subprocess.run(
            ["hmmbuild", hmm_path, sto_path],
            check=True,
            stdout=subprocess.DEVNULL,
        )

        cons_path = fam_dir / "consensus.fa"
        with open(cons_path, "w") as f:
            subprocess.run(
                ["hmmemit", "-c", hmm_path],
                check=True,
                stdout=f,
            )

        for i, t_seq in enumerate(seqs):
            domain = t_seq.extra.split()[-1]
            (pid, query_seq_name) = best_pid_by_target[domain]
            q_seq = q_sto.sequences[query_seq_name].fasta_record()

            t_seq.name += f"/{pid:.2f}/{query_seq_name}"

            query_fa_path = fam_dir / f"{i}.q.fa"
            query_hmm_path = fam_dir / f"{i}.q.hmm"
            target_fa_path = fam_dir / f"{i}.t.fa"

            with open(query_fa_path, "w") as f:
                f.write(str(q_seq))

            with open(target_fa_path, "w") as f:
                f.write(str(t_seq))

            subprocess.run(
                ["hmmbuild", query_hmm_path, target_fa_path],
                check=True,
                stdout=subprocess.DEVNULL,
            )
