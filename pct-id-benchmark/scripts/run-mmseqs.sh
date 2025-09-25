#!/usr/bin/env bash
MMSEQS=../tools/bin/mmseqs
ESL=../tools/bin/esl-seqstat

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

TMP=./tmp/mmseqs/
mkdir -p $TMP

MDB=$TMP/msaDB
QDB=$TMP/queryDB
TDB=$TMP/targetDB
ADB=$TMP/alignDB
FDB=$TMP/forwardDB

S_ARGS="$QDB $TDB $ADB $TMP --threads $THREADS -s 7.5  --max-seqs 1000"
C_ARGS="--format-mode 0"

TBL_SEQ=$RESULTS/mmseqs.seq.tbl
TIME_SEQ=$RESULTS/mmseqs.seq.time
rm -rf $TMP/*
echo "running mmseqs seq..."
(
    $MMSEQS createdb $TARGET $TDB
    $MMSEQS createdb $QUERY_FA $QDB
    $MMSEQS search $S_ARGS -e $E
    $MMSEQS convertalis $QDB $TDB $ADB $TBL_SEQ $C_ARGS
) | /usr/bin/time -p -o "$TIME_SEQ" cat >/dev/null
cat $TIME_SEQ | grep real | awk '{print $2 "s"}'
echo

TBL_PRF=$RESULTS/mmseqs.prf.tbl
TIME_PRF=$RESULTS/mmseqs.prf.time
rm -rf $TMP/*
echo "running mmseqs profile..."
(
    $MMSEQS convertmsa $QUERY_MSA $MDB --identifier-field 0
    $MMSEQS msa2profile $MDB $QDB --match-mode 1
    $MMSEQS createdb $TARGET $TDB
    $MMSEQS search $S_ARGS -e $E
    $MMSEQS convertalis $QDB $TDB $ADB $TBL_PRF $C_ARGS
) | /usr/bin/time -p -o "$TIME_PRF" cat >/dev/null
cat $TIME_PRF | grep real | awk '{print $2 "s"}'
echo

# TBL_CONS=$RESULTS/mmseqs.cons.tbl
# TIME_CONS=$RESULTS/mmseqs.cons.time
# rm -rf $TMP/*
# echo "running mmseqs consensus..."
# (
#     $MMSEQS createdb $TARGET $TDB
#     $MMSEQS createdb $QUERY_CONS_FA $QDB
#     $MMSEQS search $S_ARGS -e $E
#     $MMSEQS convertalis $QDB $TDB $ADB $TBL_CONS $C_ARGS
# ) | /usr/bin/time -p -o "$TIME_CONS" cat >/dev/null
# cat $TIME_CONS | grep real | awk '{print $2 "s"}'
# echo
