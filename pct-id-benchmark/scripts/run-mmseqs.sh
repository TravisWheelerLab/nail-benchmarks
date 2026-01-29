#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

TMP=./tmp/mmseqs/
mkdir -p $TMP

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

ARGS_SENS="-s 7.5"
ARGS_NAIL="-s 12.0 --max-seqs 2000"

QDB=$QDB_PRF
run_mmseqs "mmseqs-default.prf"
run_mmseqs "mmseqs-sens.prf" "$ARGS_SENS"
run_mmseqs "mmseqs-nail.prf" "$ARGS_NAIL"

QDB=$QDB_SEQ
run_mmseqs "mmseqs-default.seq"
run_mmseqs "mmseqs-sens.seq" "$ARGS_SENS"
run_mmseqs "mmseqs-nail.seq" "$ARGS_NAIL"

# if [ -e $P7_QDB ]; then
#     QDB=$QDB_P7
#     run_mmseqs "mmseqs-default-p7.prf"
#     run_mmseqs "mmseqs-sens-p7.prf" "$ARGS_SENS"
#     run_mmseqs "mmseqs-nail-p7.prf" "$ARGS_NAIL"
# fi
