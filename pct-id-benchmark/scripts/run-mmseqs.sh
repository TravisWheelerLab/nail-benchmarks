#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

set_default E 1e3
set_default RESULTS "${BM_DIR}/results/"

TMP=./tmp/mmseqs/
rm -rf $TMP
mkdir -p $TMP
mkdir -p $RESULTS

ANNOYING=$TMP/annoying

TDB=$TMP/targetDB
$MMSEQS createdb $TARGET $TDB > /dev/null

QDB_SEQ=$TMP/queryDB-seq
$MMSEQS createdb $QUERY_FA $QDB_SEQ > /dev/null

MDB=$TMP/msaDB
$MMSEQS convertmsa $QUERY_MSA $MDB --identifier-field 0 > /dev/null

QDB_PRF=$TMP/queryDB-prf
$MMSEQS msa2profile $MDB $QDB_PRF --match-mode 1 > /dev/null

QDB_P7=$BM_DIR/p7-queryDB/queryDB 

ADB=$TMP/alignDB

###
PREFIX_1="mmseqs-s5.7-ms2000"
PREFIX_2="mmseqs-s7.5-ms2000"
PREFIX_3="mmseqs-s10.0-ms2000"
PREFIX_4="mmseqs-s12.0-ms2000"
PREFIX_5="mmseqs-s14.0-ms2000"

ARGS_1="-s 5.7  --max-seqs 2000 -e ${E}"
ARGS_2="-s 7.5  --max-seqs 2000 -e ${E}"
ARGS_3="-s 10.0 --max-seqs 2000 -e ${E}"
ARGS_4="-s 12.0 --max-seqs 2000 -e ${E}"
ARGS_5="-s 14.0 --max-seqs 2000 -e ${E}"

QDB=$QDB_PRF
run_mmseqs "${PREFIX_1}.prf" "$ARGS_1"
run_mmseqs "${PREFIX_2}.prf" "$ARGS_2"
run_mmseqs "${PREFIX_3}.prf" "$ARGS_3"
run_mmseqs "${PREFIX_4}.prf" "$ARGS_4"
run_mmseqs "${PREFIX_5}.prf" "$ARGS_5"

QDB=$QDB_SEQ
run_mmseqs "${PREFIX_1}.seq" "$ARGS_1"
run_mmseqs "${PREFIX_2}.seq" "$ARGS_2"
run_mmseqs "${PREFIX_3}.seq" "$ARGS_3"
run_mmseqs "${PREFIX_4}.seq" "$ARGS_4"
run_mmseqs "${PREFIX_5}.seq" "$ARGS_5"


