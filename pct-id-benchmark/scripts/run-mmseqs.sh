#!/usr/bin/env bash
MMSEQS=../tools/bin/mmseqs
ESL=../tools/bin/esl-seqstat

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

TMP=./tmp-mmseqs/
mkdir -p $TMP

MDB=$TMP/msaDB
QDB=$TMP/queryDB
TDB=$TMP/targetDB
ADB=$TMP/alignDB
FDB=$TMP/forwardDB

O1=$RESULTS/mmseqs.seq.tbl
T1=$RESULTS/mmseqs.seq.time
O2=$RESULTS/mmseqs.prf.tbl
T2=$RESULTS/mmseqs.prf.time
O3=$RESULTS/mmseqs.fwd.tbl
T3=$RESULTS/mmseqs.fwd.time
O4=$RESULTS/mmseqs.fwd-max.tbl
T4=$RESULTS/mmseqs.fwd-max.time

S_ARGS="$QDB $TDB $ADB $TMP --threads $THREADS"
F_ARGS="$QDB $TDB $ADB $FDB --threads $THREADS" 
C_ARGS="--format-output target,query,tstart,tend,qstart,qend,bits,evalue"

rm -rf $TMP/*
echo "running mmseqs seq..."
(
    $MMSEQS createdb $TARGET $TDB
    $MMSEQS createdb $QUERY_FA $QDB
    $MMSEQS search $S_ARGS -e $E
    $MMSEQS convertalis $QDB $TDB $ADB $O1 $C_ARGS
) | /usr/bin/time -p -o "$T1" cat >/dev/null
cat $T1 | grep real

rm -rf $TMP/*
echo "running mmseqs profile..."
(
    $MMSEQS convertmsa $QUERY_MSA $MDB --identifier-field 0
    $MMSEQS msa2profile $MDB $QDB --match-mode 1
    $MMSEQS createdb $TARGET $TDB
    $MMSEQS search $S_ARGS -e $E
    $MMSEQS convertalis $QDB $TDB $ADB $O2 $C_ARGS
) | /usr/bin/time -p -o "$T2" cat >/dev/null
cat $T2 | grep real

rm -rf $TMP/*
N_TARGET=$($ESL "$TARGET" | grep seq | awk '{print $NF}')
PRE_E=$(echo "0.01 * $N_TARGET" | bc -l)
echo "running mmseqs fwbw..."
(
    $MMSEQS convertmsa $QUERY_MSA $MDB --identifier-field 0
    $MMSEQS msa2profile $MDB $QDB --match-mode 1
    $MMSEQS createdb $TARGET $TDB
    $MMSEQS search $S_ARGS -e $PRE_E -k 6 --k-score 80 --max-seqs 2147483647
    $MMSEQS fwbw $F_ARGS -e $E
    $MMSEQS convertalis $QDB $TDB $FDB $O3 $C_ARGS
) | /usr/bin/time -p -o "$T3" cat >/dev/null
cat $T3 | grep real

rm -rf $TMP/*
echo "running mmseqs fwbw max..."
(
    $MMSEQS convertmsa $QUERY_MSA $MDB --identifier-field 0
    $MMSEQS msa2profile $MDB $QDB --match-mode 1
    $MMSEQS createdb $TARGET $TDB
    $MMSEQS search $S_ARGS -e $E -k 6 --k-score 80 --max-seqs 2147483647
    $MMSEQS fwbw $F_ARGS -e $E
    $MMSEQS convertalis $QDB $TDB $FDB $O4 $C_ARGS
) | /usr/bin/time -p -o "$T4" cat >/dev/null
cat $T4 | grep real
