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

QUERY=$QUERY_HMM
run_nail "nail-s12.0.prf" "--allow-overwrite --mmseqs-s 12.0 -E ${E}"
run_nail "nail-s10.0.prf" "--allow-overwrite --mmseqs-s 10.0 -E ${E}"

QUERY=$QUERY_FA
run_nail "nail-s12.0.seq" "--allow-overwrite --mmseqs-s 12.0 -E ${E}"
run_nail "nail-s10.0.seq" "--allow-overwrite --mmseqs-s 10.0 -E ${E}"
