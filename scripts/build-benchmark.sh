#! /bin/sh

# * this script is expected to be called from the root of the repo
#    i.e.:
#      ./scripts/build-benchmark.sh
# *

PROFMARK_BIN=./bin/create-profmark
NAME=benchmark
N=2000000
MSA=./data/pfam.sto
FA=./data/uniprot_sprot.fasta

mkdir $NAME
$PROFMARK_BIN --single -N $N ./$NAME/$NAME $MSA $FA
hmmbuild --cpu 8 ./$NAME/$NAME.hmm ./$NAME/$NAME.msa
