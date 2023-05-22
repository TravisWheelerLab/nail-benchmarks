#! /bin/sh

NAME=benchmark
TARGET=./$NAME/$NAME.fa
QUERY_MSA=./$NAME/$NAME.msa
QUERY_HMM=./$NAME/$NAME.hmm

DIR=./output/mmoreseqs/
PREP=$DIR/prep/
SEEDS=$DIR/prep/seeds.json
OUT=$DIR/mmore.out
TSV=$DIR/mmore.tsv
SORTED=$DIR/results.sorted

EVALUE_COL=8

mkdir -p $DIR
mmoreseqs prep --skip-hmmbuild -p $PREP $QUERY_MSA $TARGET
mmoreseqs seed $PREP -s $SEEDS
mmoreseqs align -E 200 -o $TSV $QUERY_HMM $TARGET $SEEDS > $OUT

grep -v '^#' $TSV | sort -g -k$EVALUE_COL > $SORTED
