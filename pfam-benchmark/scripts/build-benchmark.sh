#! /bin/sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

PROFMARK_BIN=create-profmark
BENCHMARK_NAME=benchmark

DATA_DIR=$SCRIPT_DIR/../data/
BENCHMARK_DIR=$SCRIPT_DIR/../$BENCHMARK_NAME/

mkdir $BENCHMARK_DIR

N=2000000
MIN_TEST=10
MAX_TEST=30

MSA=$DATA_DIR/pfam.sto
FA=$DATA_DIR/uniprot_sprot.fasta

esl-sfetch --index $FA
$PROFMARK_BIN -N $N --mintest $MIN_TEST --maxtest $MAX_TEST $BENCHMARK_DIR/$BENCHMARK_NAME $MSA $FA
hmmbuild $BENCHMARK_DIR/$BENCHMARK_NAME.train.hmm $BENCHMARK_DIR/$BENCHMARK_NAME.train.msa

ln -s $DATA_DIR/long-seq $BENCHMARK_DIR/
