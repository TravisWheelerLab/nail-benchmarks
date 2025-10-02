#!/usr/bin/env bash

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

DIAMOND="$NUMA_PREFIX../tools/bin/diamond"

TMP=./tmp/diamond
mkdir -p $TMP
TARGET_DB=$TMP/target_db

TBL=$RESULTS/diamond.seq.tbl
TIME=$RESULTS/diamond.seq.time

$DIAMOND makedb --in $TARGET --db $TARGET_DB > /dev/null 2>&1

echo "running diamond seq..."
/usr/bin/time -p -o $TIME \
    $DIAMOND blastp --query $QUERY_FA \
    --db $TARGET_DB \
    --out $TBL \
    --outfmt 6 \
    --evalue $E \
    --threads $THREADS \
    > /dev/null 2>&1
cat $TIME | grep real | awk '{print $2 "s"}'
echo
