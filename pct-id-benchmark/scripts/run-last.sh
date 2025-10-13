#!/usr/bin/env bash
LASTAL="$NUMA_PREFIX../tools/bin/lastal"
LASTDB="$NUMA_PREFIX../tools/bin/lastdb"

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

TMP=./tmp/last/
mkdir -p $TMP
TARGET_DB=$TMP/target_db

$LASTDB -p $TARGET_DB $TARGET

run() {
    PREFIX=$1
    S_ARGS=$2

    TBL="$RESULTS/$PREFIX.tbl"
    TIME="$RESULTS/$PREFIX.time"

    echo "running $PREFIX | $S_ARGS"
    echo "   query: $QUERY_FA"
    echo "  target: $TARGET_DB"

    /usr/bin/time -p -o $TIME \
        $LASTAL $TARGET_DB $QUERY_FA \
        -f BlastTab \
        -P $THREADS \
        $S_ARGS \
        > $TBL

    cat $TIME | grep real | awk '{print $2 "s"}'
    echo
}

run "last.seq"
