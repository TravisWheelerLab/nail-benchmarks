#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

check_defined BM_DIR
check_defined A
check_defined B
check_defined TMP
check_defined RESULTS


TMP_NAIL="$TMP/nail"
TMP_MMSEQS="$TMP/mmseqs"
mkdir -p $TMP_NAIL
mkdir -p $TMP_MMSEQS

MGY_REV="$BM_DIR/mgy-rev/"

RESULTS_NAIL="$RESULTS/nail"
RESULTS_MMSEQS="$RESULTS/mmseqs"

mkdir -p $RESULTS_NAIL
mkdir -p $RESULTS_MMSEQS

###

QUERY=$QUERY_HMM

ANNOYING=$TMP_MMSEQS/annoying

MDB=$TMP_MMSEQS/msaDB
$MMSEQS convertmsa $QUERY_MSA $MDB --identifier-field 0 > /dev/null

QDB=$TMP_MMSEQS/queryDB-prf
$MMSEQS msa2profile $MDB $QDB --match-mode 1 > /dev/null

TDB=$TMP_MMSEQS/targetDB
ADB=$TMP_MMSEQS/alignDB

S="11.0"
MS="5000"

for i in $(seq "$A" "$B"); do
    TARGET="$MGY_REV/$i.rev.fa"

    TMP=$TMP_NAIL
    RESULTS=$RESULTS_NAIL
    run_nail "nail.${i}.rev.prf" "--allow-overwrite --mmseqs-s ${S} --mmseqs-max-seqs ${MS}"

    TMP=$TMP_MMSEQS
    RESULTS=$RESULTS_MMSEQS
    $MMSEQS createdb $TARGET $TDB > /dev/null
    run_mmseqs "mmseqs.${i}.rev.prf" "-e 100 -s ${S} --max-seqs ${MS}"
done

