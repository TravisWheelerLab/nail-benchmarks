#! /bin/sh

if [ "$#" == 0 ]; then
    echo "usage: ./run-nail.sh <benchmark-dir> [threads]"
    exit
elif [ "$#" == 1 ]; then
    THREADS=8
elif [ "$#" -ge 2 ]; then
    THREADS=$2
fi

E=1e9

BENCHMARK_DIR=$1
NAME=$(basename "$BENCHMARK_DIR")
TARGET=$BENCHMARK_DIR/$NAME.test.fa
QUERY_MSA=$BENCHMARK_DIR/$NAME.train.msa
QUERY_HMM=$BENCHMARK_DIR/$NAME.train.hmm

RESULTS_DIR=$BENCHMARK_DIR/results/nail/

PREP_TIME=$RESULTS_DIR/nail.prep.time
SEED_TIME=$RESULTS_DIR/nail.seed.time

ALIGN_DEFAULT_TIME=$RESULTS_DIR/nail.align.default.time
ALIGN_8_12_TIME=$RESULTS_DIR/nail.align.8-12.time
ALIGN_FULL_TIME=$RESULTS_DIR/nail.align.full.time

PREP=$RESULTS_DIR/prep/
PREP=$RESULTS_DIR/prep/
SEEDS=$PREP/seeds.json

OUT_DEFAULT=$RESULTS_DIR/nail.default.out
TSV_DEFAULT=$RESULTS_DIR/nail.default.tsv

OUT_8_12=$RESULTS_DIR/nail.8-12.out
TSV_8_12=$RESULTS_DIR/nail.8-12.tsv

OUT_FULL=$RESULTS_DIR/nail.full.out
TSV_FULL=$RESULTS_DIR/nail.full.tsv


rm -rf $RESULTS_DIR
mkdir -p $RESULTS_DIR
mkdir $PREP

/usr/bin/time -p -o $PREP_TIME nail prep -t $THREADS --skip-hmmbuild -p $PREP $QUERY_MSA $TARGET
/usr/bin/time -p -o $SEED_TIME nail seed -t $THREADS -q $QUERY_HMM -s $SEEDS $PREP

/usr/bin/time -p -o $ALIGN_DEFAULT_TIME nail align -t $THREADS -E $E -T $TSV_DEFAULT -O $OUT_DEFAULT $QUERY_HMM $TARGET $SEEDS
/usr/bin/time -p -o $ALIGN_8_12_TIME nail align -t $THREADS -A 8 -B 12 -E $E -T $TSV_8_12 -O $OUT_8_12 $QUERY_HMM $TARGET $SEEDS
/usr/bin/time -p -o $ALIGN_FULL_TIME nail align -t $THREADS --full-dp -E $E -T $TSV_FULL -O $OUT_FULL $QUERY_HMM $TARGET $SEEDS

