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

PREFIX_1="nail-s5.7-ms2000"
PREFIX_2="nail-s7.5-ms2000"
PREFIX_3="nail-s10.0-ms2000"
PREFIX_4="nail-s12.0-ms2000"
PREFIX_5="nail-s14.0-ms2000"

ARGS_1="--allow-overwrite --mmseqs-s 5.7  --mmseqs-max-seqs 2000 -E ${E}"
ARGS_2="--allow-overwrite --mmseqs-s 7.5  --mmseqs-max-seqs 2000 -E ${E}"
ARGS_3="--allow-overwrite --mmseqs-s 10.0 --mmseqs-max-seqs 2000 -E ${E}"
ARGS_4="--allow-overwrite --mmseqs-s 12.0 --mmseqs-max-seqs 2000 -E ${E}"
ARGS_5="--allow-overwrite --mmseqs-s 14.0 --mmseqs-max-seqs 2000 -E ${E}"

QUERY=$QUERY_HMM
run_nail "${PREFIX_1}.prf" "$ARGS_1"
run_nail "${PREFIX_2}.prf" "$ARGS_2"
run_nail "${PREFIX_3}.prf" "$ARGS_3"
run_nail "${PREFIX_4}.prf" "$ARGS_4"
run_nail "${PREFIX_5}.prf" "$ARGS_5"

QUERY=$QUERY_FA
run_nail "${PREFIX_1}.seq" "$ARGS_1"
run_nail "${PREFIX_2}.seq" "$ARGS_2"
run_nail "${PREFIX_3}.seq" "$ARGS_3"
run_nail "${PREFIX_4}.seq" "$ARGS_4"
run_nail "${PREFIX_5}.seq" "$ARGS_5"
