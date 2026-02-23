#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

set_default E 1e9
set_default THREADS_PER "2"
set_default RESULTS "${BM_DIR}/results/"

TMP=./tmp/hmmer/
rm -rf $TMP
mkdir -p $TMP
mkdir -p $RESULTS

###

QUERY=$QUERY_HMM
run_hmmsearch_split "hmmer.prf" "-E ${E}"

QUERY=$QUERY_FA
run_phmmer_split "hmmer.seq" "-E ${E}"

