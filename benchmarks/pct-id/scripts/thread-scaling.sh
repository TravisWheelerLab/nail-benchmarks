#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

THREAD_LIST=(1 2 4 8 16 32 64 96)

set_default THREAD_MAX 64
set_default E 10

check_defined BM_DIR
RESULTS="$BM_DIR/results-threads/"
mkdir -p $RESULTS

nail() {
    TMP=./tmp/nail/
    mkdir -p $TMP

    QUERY=$QUERY_HMM
    for i in "${THREAD_LIST[@]}"; do
        THREADS=$i
        run_nail "nail-t${i}-prog.prf"   "--allow-overwrite --mmseqs-s 12.0 --prog-seed"
        run_nail "nail-t${i}-ms2000.prf" "--allow-overwrite --mmseqs-s 12.0 --mmseqs-max-seqs 2000"
    done
}

hmmer() {
    TMP=./tmp/hmmer/
    QUERY=$QUERY_HMM
    for i in "${THREAD_LIST[@]}"; do
        THREADS=$i
        if (($i >= 4)); then
            THREADS_PER=2
            run_hmmsearch_split "hmmer-t${i}-spl2.prf"
        fi

        if (($i >= 8)); then
            THREADS_PER=4
            run_hmmsearch_split "hmmer-t${i}-spl4.prf"
        fi

        run_hmmsearch "hmmer-t${i}.prf"
    done
}

mmseqs() {
    TMP=./tmp/mmseqs/
    mkdir -p $TMP
    
    ANNOYING=$TMP/annoying

    MDB=$TMP/msaDB
    $MMSEQS convertmsa $QUERY_MSA $MDB --identifier-field 0 > /dev/null
    
    QDB=$TMP/queryDB-prf
    $MMSEQS msa2profile $MDB $QDB --match-mode 1 > /dev/null

    TDB=$TMP/targetDB
    $MMSEQS createdb $TARGET $TDB > /dev/null

    ADB=$TMP/alignDB
    for i in "${THREAD_LIST[@]}"; do
        THREADS=$i
        run_mmseqs "mmseqs-t${i}.prf" "-s 12.0 --max-seqs 2000"
    done
}

hmmer
nail
mmseqs
