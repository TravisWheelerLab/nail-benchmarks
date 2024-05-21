#! /bin/sh

if [ "$#" == 0 ]; then
    echo "usage: ./run-hmmsearch.sh <benchmark-dir> [threads]"
    exit
fi

if [ -n "$2" ]; then
    THREADS=$2
else
    THREADS=1
fi

E=1e9

BENCHMARK_DIR=$1
NAME=$(basename "$BENCHMARK_DIR")
TARGET=$BENCHMARK_DIR/$NAME.test.fa
QUERY=$BENCHMARK_DIR/$NAME.train.hmm

RESULTS_DIR=$BENCHMARK_DIR/results/hmmer/

TIME=$RESULTS_DIR/hmmer.time

OUT=$RESULTS_DIR/hmmer.out
TBL=$RESULTS_DIR/hmmer.tbl
DOM=$RESULTS_DIR/hmmer.domtbl

if [ -d "$RESULTS_DIR" ]; then
    echo "results directory already exists"
    exit
fi

mkdir -p $RESULTS_DIR

echo "running hmmsearch..."
/usr/bin/time -p -o $TIME \
    hmmsearch \
    --cpu $THREADS \
    -E $E \
    -o $OUT \
    --notextw \
    --domtblout $DOM \
    --tblout $TBL \
    $QUERY $TARGET

awk '/real/ {print "time:", $2}' $TIME
