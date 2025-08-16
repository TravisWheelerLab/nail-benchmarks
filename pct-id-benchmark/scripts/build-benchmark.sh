#! /bin/sh

START_DIR="$(pwd)"
SCRIPT_DIR="$(cd -- "$(dirname -- "$0")" && pwd)"
PROFMARK_BIN=create-profmark
NAME=benchmark
QUERY_HMM=$NAME.train.hmm
QUERY_MSA=$NAME.train.msa
TARGET_FA=$NAME.test.fa

TRAIN_TEST_ID=0.5
MIN_TEST=10
MAX_TEST=30

if [ "$#" == 0 ]; then
    echo "usage: ./build-benchmark.sh <source.sto> <random.fasta> <n-decoys> <dir/>"
    exit 1
elif [ "$#" == 4 ]; then
    MSA=$(realpath "$1")
    FA=$(realpath "$2")
    N_DECOYS=$3
    DIR="$4"
fi

if [ -z "$DIR" ]; then
    echo "error: <dir/> is: '$DIR'"
    exit 1
fi

PM_DIR=$DIR/profmark/
mkdir -p $PM_DIR

esl-sfetch --index $FA >> /dev/null

$PROFMARK_BIN \
    -N $N_DECOYS \
    -1 $TRAIN_TEST_ID \
    --mintest $MIN_TEST \
    --maxtest $MAX_TEST \
    $PM_DIR/$NAME $MSA $FA

status=$?
if [ $status -ne 0 ]; then
    echo "\033[31mprofmark failed\033[0m"
    exit 1
fi

cd $PM_DIR && \
    ln -sf $MSA source.sto && \
    ln -sf $QUERY_MSA query.sto && \
    ln -sf $TARGET_FA target.fa &&
    
cd $START_DIR

python3 $SCRIPT_DIR/sample.py $DIR
