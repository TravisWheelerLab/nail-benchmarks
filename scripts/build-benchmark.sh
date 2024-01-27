#! /bin/sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

PROFMARK_BIN=../../bin/profmark-3.4
N=2000000
MIN_TEST=10
MAX_TEST=30

MSA=../data/pfam.sto
FA=../data/uniprot_sprot.fasta

NAME=benchmark
$PROFMARK_BIN -N $N --mintest $MIN_TEST --maxtest $MAX_TEST ./$NAME/$NAME  $MSA $FA
hmmbuild ./$NAME/$NAME.train.hmm ./$NAME/$NAME.train.msa
