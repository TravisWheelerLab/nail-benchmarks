#! /bin/sh

if [ "$#" == 0 ]; then
    echo "usage: ./run-nail.sh <benchmark-dir> [threads]"
    exit
elif [ "$#" == 1 ]; then
    THREADS=8
elif [ "$#" -ge 2 ]; then
    THREADS=$2
fi

HMMSEARCH=hmmsearch

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

/usr/bin/time -p -o $TIME $HMMSEARCH --cpu $THREADS -E 200 -o $OUT --domtblout $DOM --tblout $TBL $QUERY $TARGET
