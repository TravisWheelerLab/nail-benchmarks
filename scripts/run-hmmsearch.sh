#! /bin/sh

if [ "$#" == 0 ]; then
    echo "usage: ./run-nail.sh <benchmark-dir>"
    exit
fi

E=1e9

BENCHMARK_DIR=$1
NAME=$(basename "$BENCHMARK_DIR")
TARGET=$BENCHMARK_DIR/$NAME.test.fa
QUERY=$BENCHMARK_DIR/$NAME.train.hmm

RESULTS_DIR=$BENCHMARK_DIR/results/hmmer/

TIME_1=$RESULTS_DIR/hmmer.1.time
TIME_8=$RESULTS_DIR/hmmer.8.time

OUT=$RESULTS_DIR/hmmer.out
TBL=$RESULTS_DIR/hmmer.tbl
DOM=$RESULTS_DIR/hmmer.domtbl

if [ -d "$RESULTS_DIR" ]; then
    echo "results directory already exists"
    exit
fi

mkdir -p $RESULTS_DIR

/usr/bin/time -p -o $TIME_1 hmmsearch --cpu 1 -E $E -o $OUT --domtblout $DOM --tblout $TBL $QUERY $TARGET
/usr/bin/time -p -o $TIME_8 hmmsearch --cpu 8 -E $E -o $OUT --domtblout $DOM --tblout $TBL $QUERY $TARGET
