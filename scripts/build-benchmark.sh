#! /bin/sh

# * this script is expected to be called from the root of the repo
#    i.e.:
#      ./scripts/build-benchmark.sh
# *

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

PROFMARK_BIN=./bin/create-profmark
NAME=benchmark
N=2000000
MSA=./data/pfam.sto
FA=./data/uniprot_sprot.fasta

mkdir $NAME
$PROFMARK_BIN -N $N ./$NAME/$NAME $MSA $FA
hmmbuild --cpu 24 ./$NAME/$NAME.train.hmm ./$NAME/$NAME.msa
mv ./$NAME/$NAME.test.fa ./$NAME/$NAME.fa
