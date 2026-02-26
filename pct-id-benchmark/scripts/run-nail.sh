#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

set_default E 1e9
set_default RESULTS "${BM_DIR}/results/"

TMP=./tmp/nail/
rm -rf $TMP
mkdir -p $TMP
mkdir -p $RESULTS

###

ARGS_1="--allow-overwrite --mmseqs-s 10.0 --prog-seed            -E ${E}"
ARGS_2="--allow-overwrite --mmseqs-s 10.0 --mmseqs-max-seqs 2000 -E ${E}"

QUERY=$QUERY_HMM
run_nail "nail-prog.prf" "$ARGS_1"
run_nail "nail-ms2000.prf" "$ARGS_2"

#QUERY=$QUERY_FA
#run_nail "nail-s12.0.seq" "$ARGS_1"
#run_nail "nail-s10.0.seq" "$ARGS_2"
