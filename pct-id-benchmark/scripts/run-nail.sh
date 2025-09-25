#!/usr/bin/env bash
NAIL=../tools/bin/nail

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

TMP=./tmp/nail/
mkdir -p $TMP

STATS=$RESULTS/nail.stats

S_ARGS="-s -t $THREADS -E $E --tmp-dir $TMP --double-seed"

TBL_SEQ=$RESULTS/nail.seq.tbl
TIME_SEQ=$RESULTS/nail.seq.time
echo "running nail seq..."
/usr/bin/time -p -o $TIME_SEQ \
    $NAIL search $S_ARGS \
    --tbl-out $TBL_SEQ \
    $QUERY_FA $TARGET >> $STATS
cat $TIME_SEQ | grep real | awk '{print $2 "s"}'
echo

TBL_PRF=$RESULTS/nail.prf.tbl
TIME_PRF=$RESULTS/nail.prf.time
echo "running nail hmm..."
/usr/bin/time -p -o $TIME_PRF \
    $NAIL search $S_ARGS \
    --tbl-out $TBL_PRF \
    $QUERY_HMM $TARGET >> $STATS
cat $TIME_PRF | grep real | awk '{print $2 "s"}'
echo

# TBL_CONS=$RESULTS/nail.cons.tbl
# TIME_CONS=$RESULTS/nail.cons.time
# echo "running nail consensus..."
# /usr/bin/time -p -o $TIME_CONS \
#     $NAIL search $S_ARGS \
#     --tbl-out $TBL_CONS \
#     $QUERY_CONS_FA $TARGET >> $STATS
# cat $TIME_CONS | grep real | awk '{print $2 "s"}'
# echo

