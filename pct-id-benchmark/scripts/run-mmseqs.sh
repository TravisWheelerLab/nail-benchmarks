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

TBL_1=$RESULTS/mmseqs.seq.tbl
TIME_1=$RESULTS/mmseqs.seq.time
TBL_2=$RESULTS/mmseqs.prf.tbl
TIME_2=$RESULTS/mmseqs.prf.time

S_ARGS="$QDB $TDB $ADB $TMP --threads $THREADS"
C_ARGS="--format-mode 0"

rm -rf $TMP/*
echo "running mmseqs seq..."
(
    $MMSEQS createdb $TARGET $TDB
    $MMSEQS createdb $QUERY_FA $QDB
    $MMSEQS search $S_ARGS -e $E
    $MMSEQS convertalis $QDB $TDB $ADB $TBL_1 $C_ARGS
) | /usr/bin/time -p -o "$TIME_1" cat >/dev/null
cat $TIME_1 | grep real | awk '{print $2 "s"}'
echo

rm -rf $TMP/*
echo "running mmseqs profile..."
(
    $MMSEQS convertmsa $QUERY_MSA $MDB --identifier-field 0
    $MMSEQS msa2profile $MDB $QDB --match-mode 1
    $MMSEQS createdb $TARGET $TDB
    $MMSEQS search $S_ARGS -e $E
    $MMSEQS convertalis $QDB $TDB $ADB $TBL_2 $C_ARGS
) | /usr/bin/time -p -o "$TIME_2" cat >/dev/null
cat $TIME_2 | grep real | awk '{print $2 "s"}'
echo
