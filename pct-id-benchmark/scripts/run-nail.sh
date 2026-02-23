#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

TMP=./tmp/nail/
mkdir -p $TMP

###

QUERY=$QUERY_HMM
run_nail "nail-s12.0.prf"       "--allow-overwrite --mmseqs-s 12.0"
run_nail "nail-s10.0.prf"       "--allow-overwrite --mmseqs-s 10.0"

QUERY=$QUERY_FA
run_nail "nail-s12.0.seq" "--mmseqs-s 12.0"
run_nail "nail-s10.0.seq" "--mmseqs-s 10.0"
