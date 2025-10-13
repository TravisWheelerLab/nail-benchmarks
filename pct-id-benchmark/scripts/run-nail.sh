#!/usr/bin/env bash

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

NAIL="$NUMA_PREFIX../tools/bin/nail"

TMP=./tmp/nail/
mkdir -p $TMP

SUMMARY=$RESULTS/nail.summary

run() {
    PREFIX=$1
    S_ARGS=$2

    TBL="$RESULTS/$PREFIX.tbl"
    TIME="$RESULTS/$PREFIX.time"
    SEEDS="$RESULTS/$PREFIX.seeds"
    STATS="$RESULTS/$PREFIX.stats"

    echo "running $PREFIX | $S_ARGS"
    echo "   query: $QUERY"
    echo "  target: $TARGET"

    /usr/bin/time -p -o $TIME \
        $NAIL search \
        -s \
        -t $THREADS \
        --tmp-dir $TMP \
        --stats-results-path $STATS \
        --tbl-out $TBL \
        -E $E \
        $S_ARGS \
        $QUERY $TARGET >> $SUMMARY

    mv $TMP/align_a.tsv $SEEDS

    cat $TIME | grep real | awk '{print $2 "s"}'
    echo
}


QUERY=$QUERY_HMM
run "nail.prf" "--mmseqs-k-score 60 --mmseqs-max-seqs 2000 -C 0.01"

# copy the p7hmm-derived mmseqs2 
# profile DB to the benchmark dir
mkdir -p $DIR/p7-queryDB
cp $TMP/queryDB* $DIR/p7-queryDB

QUERY=$QUERY_FA
run "nail.seq" "--mmseqs-k-score 60 --mmseqs-max-seqs 2000 -C 0.01"
