#! /bin/sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

PROFMARK_BIN=create-profmark
BENCHMARK_NAME=benchmark

N=2000000
MIN_TEST=10
MAX_TEST=30

MSA=../data/pfam.sto
FA=../data/uniprot_sprot.fasta

$PROFMARK_BIN -N $N --mintest $MIN_TEST --maxtest $MAX_TEST ./$BENCHMARK_NAME/$BENCHMARK_NAME  $MSA $FA
hmmbuild ./$BENCHMARK_NAME/$BENCHMARK_NAME.train.hmm ./$BENCHMARK_NAME/$BENCHMARK_NAME.train.msa
