#!/usr/bin/env bash

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

NAIL="$NUMA_PREFIX../tools/bin/nail"

TMP=./tmp/nail/
mkdir -p $TMP

STATS=$RESULTS/nail.stats

S_ARGS="-s -t $THREADS -E $E --tmp-dir $TMP --mmseqs-k-score 60"

PREFIX="nail.seq"
TBL="$RESULTS/$PREFIX.tbl"
TIME="$RESULTS/$PREFIX.time"
echo "running $PREFIX..."
/usr/bin/time -p -o $TIME \
    $NAIL search $S_ARGS \
    --tbl-out $TBL \
    $QUERY_FA $TARGET >> $STATS
cat $TIME | grep real | awk '{print $2 "s"}'
echo

PREFIX="nail.prf"
TBL="$RESULTS/$PREFIX.tbl"
TIME="$RESULTS/$PREFIX.time"
echo "running $PREFIX..."
/usr/bin/time -p -o $TIME \
    $NAIL search $S_ARGS \
    --tbl-out $TBL \
    $QUERY_HMM $TARGET >> $STATS
cat $TIME | grep real | awk '{print $2 "s"}'
echo

# copy the p7hmm-derived mmseqs2 
# profile DB to the benchmark dir
mkdir -p $DIR/p7-queryDB
cp $TMP/queryDB* $DIR/p7-queryDB

