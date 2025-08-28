#!/usr/bin/env bash
LASTAL=../tools/bin/lastal
LASTDB=../tools/bin/lastdb

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

TMP=./tmp/last/
mkdir -p $TMP
TARGET_DB=$TMP/target_db

TBL_1=$RESULTS/last.seq.tbl
TIME_1=$RESULTS/last.seq.time

$LASTDB -p $TARGET_DB $TARGET

echo "running last..."
/usr/bin/time -p -o $TIME_1 \
    $LASTAL $TARGET_DB $QUERY_FA \
    -f BlastTab \
    -P $THREADS \
    > $TBL_1
cat $TIME_1 | grep real | awk '{print $2 "s"}'
echo
