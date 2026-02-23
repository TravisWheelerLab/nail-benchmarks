#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

set_default E 1e9
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

ARGS_DEFAULT="-e ${E}"
ARGS_SENS="-s 7.5 -e ${E}"

QDB=$QDB_PRF
run_mmseqs "mmseqs-default.prf" "$ARGS_DEFAULT"
run_mmseqs "mmseqs-sens.prf" "$ARGS_SENS"

QDB=$QDB_SEQ
run_mmseqs "mmseqs-default.seq" "$ARGS_DEFAULT"
run_mmseqs "mmseqs-sens.seq" "$ARGS_SENS"


