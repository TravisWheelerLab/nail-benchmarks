#!/usr/bin/env bash
NAIL=../tools/bin/nail

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

echo $THREADS

TMP=./tmp-nail/
mkdir -p $TMP

STATS=$RESULTS/nail.stats
O1=$RESULTS/nail.seq.tbl
T1=$RESULTS/nail.seq.time
O2=$RESULTS/nail.hmm.tbl
T2=$RESULTS/nail.hmm.time
O3=$RESULTS/nail.hmm-double.tbl
T3=$RESULTS/nail.hmm-double.time

S_ARGS="-s -t $THREADS -E $E --tmp-dir $TMP"

echo "running nail seq..."
/usr/bin/time -p -o $T1 \
    $NAIL search $S_ARGS \
    --tbl-out $O1 \
    $QUERY_FA $TARGET >> $STATS
cat $T1 | grep real

echo "running nail hmm..."
/usr/bin/time -p -o $T2 \
    $NAIL search $S_ARGS \
    --tbl-out $O2 \
    $QUERY_HMM $TARGET >> $STATS
cat $T2 | grep real

# echo "running nail hmm double..."
# /usr/bin/time -p -o $T3 \
#     nail search $S_ARGS \
#     --double-seed \
#     --tbl-out $O3 \
#     $QUERY_HMM $TARGET >> $STATS
# cat $T3 | grep real
