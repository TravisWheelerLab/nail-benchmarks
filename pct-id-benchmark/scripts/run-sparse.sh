#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

TMP=./tmp/nail/
mkdir -p $TMP

###

E=1e9
set_default RESULTS "$BM_DIR/results-sparse/"

QUERY=$QUERY_HMM
run_nail "nail-s12.0-cells.prf" "--allow-overwrite --mmseqs-s 12.0 --mmseqs-max-seqs 200 -S ${E} -C ${E} -F ${E} -E ${E} --f32-p 5"
run_nail "nail-s12.0-full.prf"  "--allow-overwrite --mmseqs-s 12.0 --mmseqs-max-seqs 200 -S ${E} -C ${E} -F ${E} -E ${E} --full-dp"

