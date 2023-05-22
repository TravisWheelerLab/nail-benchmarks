#! /bin/sh

NAME=benchmark

DIR=./output/mmoreseqs/
PREP=$DIR/prep/
OUT=$DIR/mmore.out
TSV=$DIR/mmore.tsv
SORTED=$DIR/results.sorted

EVALUE_COL=8

mkdir -p $DIR
mmoreseqs search -E 200 -p $PREP -o $TSV ./$NAME/$NAME.msa ./$NAME/$NAME.fa > $OUT

grep -v '^#' $TSV | sort -g -k$EVALUE_COL > $SORTED
