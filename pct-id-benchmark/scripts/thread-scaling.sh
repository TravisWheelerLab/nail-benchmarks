#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

set_default THREAD_MAX 32
set_default E 10

check_defined BM_DIR
RESULTS="$BM_DIR/results-threads/"
mkdir -p $RESULTS

nail() {
    TMP=./tmp/nail/
    mkdir -p $TMP

    QUERY=$QUERY_HMM
    for ((i=1; i<=$THREAD_MAX; i++)); do
        THREADS=$i
        run_nail "nail-t${i}.prf" "--mmseqs-s 12.0 --mmseqs-max-seqs 2000 -C 0.01"
    done
}

hmmer() {
    QUERY=$QUERY_HMM
    for ((i=1; i<=$THREAD_MAX; i++)); do
        THREADS=$i
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
    
    for ((i=1; i<=$THREAD_MAX; i++)); do
        THREADS=$i
        run_mmseqs "mmseqs-t${i}.prf" "-s 12.0 --max-seqs 2000"
    done
}

mmseqs
nail
hmmer
