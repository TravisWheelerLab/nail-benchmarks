#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

TMP=./tmp/hmmer/
mkdir -p $TMP

###

QUERY=$QUERY_HMM
run_hmmsearch_split "hmmer.prf"

QUERY=$QUERY_FA
run_phmmer_split "hmmer.seq"

