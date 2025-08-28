#!/usr/bin/env bash
DIAMOND=../tools/bin/diamond

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

TMP=./tmp/diamond
mkdir -p $TMP
TARGET_DB=$TMP/target_db

TBL_1=$RESULTS/diamond.seq.tbl
TIME_1=$RESULTS/diamond.seq.time

$DIAMOND makedb --in $TARGET --db $TARGET_DB > /dev/null 2>&1

echo "running diamond..."
/usr/bin/time -p -o $TIME_1 \
    $DIAMOND blastp --query $QUERY_FA \
    --db $TARGET_DB \
    --out $TBL_1 \
    --outfmt 6 \
    --evalue $E \
    --threads $THREADS \
    > /dev/null 2>&1
cat $TIME_1 | grep real | awk '{print $2 "s"}'
echo
