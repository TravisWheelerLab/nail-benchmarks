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
P7_QDB=$DIR/p7-queryDB/queryDB 

S_ARGS="$QDB $TDB $ADB $TMP --threads $THREADS -s 7.5  --max-seqs 1000"
C_ARGS="--format-mode 0"

TBL=$RESULTS/mmseqs.seq.tbl
TIME=$RESULTS/mmseqs.seq.time

rm -rf $TMP/*
echo "running mmseqs seq..."
$MMSEQS createdb $TARGET $TDB > /dev/null
$MMSEQS createdb $QUERY_FA $QDB > /dev/null
(
    $MMSEQS search $S_ARGS -e $E
    $MMSEQS convertalis $QDB $TDB $ADB $TBL $C_ARGS
) | /usr/bin/time -p -o "$TIME" cat >/dev/null
cat $TIME | grep real | awk '{print $2 "s"}'
echo

TBL=$RESULTS/mmseqs.prf.tbl
TIME=$RESULTS/mmseqs.prf.time
rm -rf $TMP/*
echo "running mmseqs profile..."
$MMSEQS createdb $TARGET $TDB > /dev/null
$MMSEQS convertmsa $QUERY_MSA $MDB --identifier-field 0 > /dev/null
$MMSEQS msa2profile $MDB $QDB --match-mode 1 > /dev/null
(
    $MMSEQS search $S_ARGS -e $E
    $MMSEQS convertalis $QDB $TDB $ADB $TBL $C_ARGS
) | /usr/bin/time -p -o "$TIME" cat >/dev/null
cat $TIME | grep real | awk '{print $2 "s"}'
echo

if [ -e $P7_QDB ]; then
    TBL=$RESULTS/mmseqs-p7.prf.tbl
    TIME=$RESULTS/mmseqs-p7.prf.time
    rm -rf $TMP/*
    echo "running mmseqs p7 profile..."
    $MMSEQS createdb $TARGET $TDB > /dev/null
    (
        $MMSEQS search $P7_QDB $TDB $ADB $TMP --threads $THREADS -s 7.5  --max-seqs 1000
        $MMSEQS convertalis $P7_QDB $TDB $ADB $TBL $C_ARGS
    ) | /usr/bin/time -p -o "$TIME" cat >/dev/null
    cat $TIME | grep real | awk '{print $2 "s"}'
    echo
fi

# TBL=$RESULTS/mmseqs.cons.tbl
# TIME=$RESULTS/mmseqs.cons.time
# rm -rf $TMP/*
# echo "running mmseqs consensus..."
# $MMSEQS createdb $QUERY_CONS_FA $QDB
# (
#     $MMSEQS search $S_ARGS -e $E
#     $MMSEQS convertalis $QDB $TDB $ADB $TBL $C_ARGS
# ) | /usr/bin/time -p -o "$TIME" cat >/dev/null
# cat $TIME | grep real | awk '{print $2 "s"}'
# echo
