#!/usr/bin/env bash

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

MMSEQS="$NUMA_PREFIX../tools/bin/mmseqs"

TMP=./tmp/mmseqs/
mkdir -p $TMP

ANNOYING=$TMP/annoying

TDB=$TMP/targetDB
$MMSEQS createdb $TARGET $TDB > /dev/null

QDB_SEQ=$TMP/queryDB-seq
$MMSEQS createdb $QUERY_FA $QDB_SEQ > /dev/null

MDB=$TMP/msaDB
$MMSEQS convertmsa $QUERY_MSA $MDB --identifier-field 0 > /dev/null

QDB_PRF=$TMP/queryDB-prf
$MMSEQS msa2profile $MDB $QDB_PRF --match-mode 1 > /dev/null

# QDB_PRF_NC=$TMP/queryDB-prf-no-comp
# $MMSEQS msa2profile $MDB $QDB_PRF_NC --match-mode 1 --comp-bias-corr 0 > /dev/null

QDB_P7=$DIR/p7-queryDB/queryDB 

ADB=$TMP/alignDB

run() {
    PREFIX=$1
    S_ARGS=$2

    [ -e $ANNOYING ] && rm -rf $ANNOYING
    [ -e $ADB.1 ] && rm -f $ADB.*

    TBL="$RESULTS/$PREFIX.tbl"
    TIME="$RESULTS/$PREFIX.time"
    echo "running $PREFIX | $S_ARGS"
    echo "   query: $QDB"
    echo "  target: $TDB"
    /usr/bin/time -p -o $TIME \
        $MMSEQS search \
        $QDB $TDB $ADB $ANNOYING \
        --threads $THREADS \
        -e $E \
        $S_ARGS > /dev/null

    $MMSEQS convertalis $QDB $TDB $ADB $TBL --format-mode 0 > /dev/null

    cat $TIME | grep real | awk '{print $2 "s"}'
    echo
}

ARGS_SENS="-s 7.5"
ARGS_NAIL="--k-score 60 --max-seqs 2000"

QDB=$QDB_PRF
run "mmseqs-default.prf"
run "mmseqs-sens.prf" "$ARGS_SENS"
run "mmseqs-nail.prf" "$ARGS_NAIL"


QDB=$QDB_SEQ
run "mmseqs-default.seq"
run "mmseqs-sens.seq" "$ARGS_SENS"
run "mmseqs-nail.seq" "$ARGS_NAIL"


if [ -e $P7_QDB ]; then
    QDB=$QDB_P7
    run "mmseqs-default-p7.prf"
    run "mmseqs-sens-p7.prf" "$ARGS_SENS"
    run "mmseqs-nail-p7.prf" "$ARGS_NAIL"
fi
