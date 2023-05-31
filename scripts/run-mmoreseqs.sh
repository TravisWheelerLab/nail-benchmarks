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
/usr/bin/time -p -o $PREP_TIME mmoreseqs prep -t $THREADS --skip-hmmbuild -p $PREP $QUERY_MSA $TARGET
/usr/bin/time -p -o $SEED_TIME mmoreseqs seed -t $THREADS -q $QUERY_HMM -s $SEEDS $PREP
/usr/bin/time -p -o $ALIGN_TIME mmoreseqs align -t $THREADS -E 200 -T $TSV $QUERY_HMM $TARGET $SEEDS

sort -g -k$EVALUE_COL $TSV > $SORTED
