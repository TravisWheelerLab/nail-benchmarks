#! /bin/sh

if [ "$#" == 0 ]; then
    echo "usage: ./run-all.sh <benchmark-dir> [threads]"
    exit
fi

if [ -n "$2" ]; then
    THREADS=$2
else
    THREADS=1
fi

BENCHMARK_DIR=$1
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

$SCRIPT_DIR/run-nail.sh $BENCHMARK_DIR $THREADS
$SCRIPT_DIR/run-mmseqs.sh $BENCHMARK_DIR $THREADS
$SCRIPT_DIR/run-hmmsearch.sh $BENCHMARK_DIR $THREADS
