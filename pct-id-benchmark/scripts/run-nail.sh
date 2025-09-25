#!/usr/bin/env bash
NAIL=../tools/bin/nail

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

TMP=./tmp/nail/
mkdir -p $TMP

STATS=$RESULTS/nail.stats

S_ARGS="-s -t $THREADS -E $E --tmp-dir $TMP"

TBL=$RESULTS/nail.seq.tbl
TIME=$RESULTS/nail.seq.time
echo "running nail seq..."
/usr/bin/time -p -o $TIME \
    $NAIL search $S_ARGS \
    --tbl-out $TBL \
    $QUERY_FA $TARGET >> $STATS
cat $TIME | grep real | awk '{print $2 "s"}'
echo

TBL=$RESULTS/nail-dbl.seq.tbl
TIME=$RESULTS/nail-dbl.seq.time
echo "running nail seq..."
/usr/bin/time -p -o $TIME \
    $NAIL search $S_ARGS --double-seed \
    --tbl-out $TBL \
    $QUERY_FA $TARGET >> $STATS
cat $TIME | grep real | awk '{print $2 "s"}'
echo

TBL=$RESULTS/nail.prf.tbl
TIME=$RESULTS/nail.prf.time
echo "running nail hmm..."
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

TBL=$RESULTS/nail-dbl.prf.tbl
TIME=$RESULTS/nail-dbl.prf.time
echo "running nail hmm double seed..."
/usr/bin/time -p -o $TIME \
    $NAIL search $S_ARGS --double-seed \
    --tbl-out $TBL \
    $QUERY_HMM $TARGET >> $STATS
cat $TIME | grep real | awk '{print $2 "s"}'
echo

# TBL=$RESULTS/nail.cons.tbl
# TIME=$RESULTS/nail.cons.time
# echo "running nail consensus..."
# /usr/bin/time -p -o $TIME \
#     $NAIL search $S_ARGS \
#     --tbl-out $TBL \
#     $QUERY_CONS_FA $TARGET >> $STATS
# cat $TIME | grep real | awk '{print $2 "s"}'
# echo

