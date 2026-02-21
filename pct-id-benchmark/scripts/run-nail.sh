#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

TMP=./tmp/nail/
mkdir -p $TMP

###


QUERY=$QUERY_HMM
run_nail "nail-s12.0.prf"       "--allow-overwrite --mmseqs-s 12.0"
# run_nail "nail-s10.0.prf"       "--allow-overwrite --mmseqs-s 10.0"
# run_nail "nail-s7.5.prf"        "--allow-overwrite --mmseqs-s 7.5"

# copy the p7hmm-derived mmseqs2 
# profile DB to the benchmark dir
mkdir -p $BM_DIR/p7-queryDB
cp $TMP/queryDB* $BM_DIR/p7-queryDB

QUERY=$QUERY_FA
run_nail "nail-s12.0.seq" "--mmseqs-s 12.0"
# run_nail "nail-s10.0.seq" "--mmseqs-s 10.0"
# run_nail "nail-s7.5.seq"  "--mmseqs-s 7.5"
