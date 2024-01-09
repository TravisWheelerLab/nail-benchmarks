#! /bin/sh

if [ "$#" == 0 ]; then
    THREADS=8
elif [ "$#" -ge 1 ]; then
    THREADS=$1
fi

PWD=$(pwd)
NAME=benchmark
TARGET=./$NAME/$NAME.fa
QUERY_MSA=./$NAME/$NAME.msa
QUERY_HMM=$PWD/$NAME/$NAME.hmm

DIR=./output/nail/
PREP_TIME=$DIR/nail.prep.time
SEED_TIME=$DIR/nail.seed.time
ALIGN_TIME=$DIR/nail.align.time

PREP=$DIR/prep/
SEEDS=$DIR/prep/seeds.json
OUT=$DIR/nail.out
TSV=$DIR/nail.tsv
SORTED=$DIR/results.sorted

EVALUE_COL=8

mkdir -p $DIR
/usr/bin/time -p -o $PREP_TIME nail prep -t $THREADS --skip-hmmbuild -p $PREP $QUERY_MSA $TARGET
/usr/bin/time -p -o $SEED_TIME nail seed -t $THREADS -q $QUERY_HMM -s $SEEDS $PREP
/usr/bin/time -p -o $ALIGN_TIME nail align -t $THREADS -E 200 -T $TSV -O $OUT $QUERY_HMM $TARGET $SEEDS

sort -g -k$EVALUE_COL $TSV > $SORTED
