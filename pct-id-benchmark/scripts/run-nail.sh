#!/usr/bin/env bash
NAIL=../tools/bin/nail

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

TMP=./tmp/nail/
mkdir -p $TMP

STATS=$RESULTS/nail.stats
TBL_1=$RESULTS/nail.pair.tbl
TIME_1=$RESULTS/nail.pair.time
TBL_2=$RESULTS/nail.cons.tbl
TIME_2=$RESULTS/nail.cons.time
TBL_3=$RESULTS/nail.prf.tbl
TIME_3=$RESULTS/nail.prf.time

S_ARGS="-s -t $THREADS -E $E --tmp-dir $TMP --double-seed"

echo "running nail pairwise..."
/usr/bin/time -p -o $TIME_1 \
    $NAIL search $S_ARGS \
    --tbl-out $TBL_1 \
    $QUERY_FA $TARGET >> $STATS
cat $TIME_1 | grep real | awk '{print $2 "s"}'
echo

echo "running nail consensus..."
/usr/bin/time -p -o $TIME_2 \
    $NAIL search $S_ARGS \
    --tbl-out $TBL_2 \
    $QUERY_CONS_FA $TARGET >> $STATS
cat $TIME_2 | grep real | awk '{print $2 "s"}'
echo

echo "running nail hmm..."
/usr/bin/time -p -o $TIME_3 \
    $NAIL search $S_ARGS \
    --tbl-out $TBL_3 \
    $QUERY_HMM $TARGET >> $STATS
cat $TIME_3 | grep real | awk '{print $2 "s"}'
echo
