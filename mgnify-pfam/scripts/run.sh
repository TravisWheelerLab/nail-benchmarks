#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

set_default THREADS 8
set_default E 10

check_defined BM_DIR
RESULTS="$BM_DIR/results/"
mkdir -p $RESULTS

nail() {
    TMP=./tmp/nail/
    QUERY=$QUERY_HMM
    run_nail "nail.prf" "--mmseqs-s 12.0 --mmseqs-max-seqs 2000 -C 0.01"
}

hmmer() {
    QUERY=$QUERY_HMM
    run_hmmsearch "hmmer.prf"
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
    
    run_mmseqs "mmseqs.prf" "-s 12.0 --max-seqs 2000"
}

hmmer
nail
mmseqs
