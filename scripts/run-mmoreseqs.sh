#! /bin/sh

PWD=$(pwd)
NAME=benchmark
TARGET=./$NAME/$NAME.fa
QUERY_MSA=./$NAME/$NAME.msa
QUERY_HMM=$PWD/$NAME/$NAME.hmm

DIR=./output/mmoreseqs/
PREP_TIME=$DIR/mmore.prep.time
SEED_TIME=$DIR/mmore.seed.time
ALIGN_TIME=$DIR/mmore.align.time

PREP=$DIR/prep/
SEEDS=$DIR/prep/seeds.json
OUT=$DIR/mmore.out
TSV=$DIR/mmore.tsv
SORTED=$DIR/results.sorted

EVALUE_COL=8

mkdir -p $DIR
/usr/bin/time -o $PREP_TIME mmoreseqs prep --skip-hmmbuild -p $PREP $QUERY_MSA $TARGET
/usr/bin/time -o $SEED_TIME mmoreseqs seed --query-hmm $QUERY_HMM -s $SEEDS $PREP
/usr/bin/time -o $ALIGN_TIME mmoreseqs align -E 200 -o $TSV $QUERY_HMM $TARGET $SEEDS > $OUT

grep -v '^#' $TSV | sort -g -k$EVALUE_COL > $SORTED
