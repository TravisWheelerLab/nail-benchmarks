#! /bin/sh

NAME=benchmark
TARGET=./$NAME/$NAME.fa
QUERY=./$NAME/$NAME.msa

DIR=./output/mmoreseqs/
PREP=$DIR/prep/
OUT=$DIR/mmore.out
TSV=$DIR/mmore.tsv
SORTED=$DIR/results.sorted

EVALUE_COL=8

mkdir -p $DIR
mmoreseqs search -E 200 -p $PREP -o $TSV $QUERY $TARGET > $OUT

grep -v '^#' $TSV | sort -g -k$EVALUE_COL > $SORTED
