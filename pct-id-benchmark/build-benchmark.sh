#! /bin/sh

PROFMARK_BIN=create-profmark
BENCHMARK_NAME=benchmark


if [ "$#" == 0 ]; then
    echo "usage: ./build-benchmark.sh <queries.sto> <random.fasta> <n-decoys> <dir/>"
    exit
elif [ "$#" == 4 ]; then
    MSA=$1
    FA=$2
    N=$3
    DIR=$4
fi

mkdir $DIR

QUERY_HMM=$BENCHMARK_NAME.train.hmm
QUERY_MSA=$BENCHMARK_NAME.train.msa
TARGET_FA=$BENCHMARK_NAME.test.fa

TRAIN_TEST_ID=0.5
MIN_TEST=10
MAX_TEST=30

esl-sfetch --index $FA

$PROFMARK_BIN -N $N \
    -1 $TRAIN_TEST_ID \
    --mintest $MIN_TEST \
    --maxtest $MAX_TEST \
    $DIR/$BENCHMARK_NAME $MSA $FA

cd $DIR && \
    hmmbuild -- cpu 8 \
    $QUERY_HMM $QUERY_MSA && \
    ln -s $QUERY_MSA query.sto && \
    ln -s $QUERY_HMM query.hmm && \
    ln -s $TARGET_FA target.fa
