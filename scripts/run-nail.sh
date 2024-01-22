#! /bin/sh

if [ "$#" == 0 ]; then
    echo "usage: ./run-nail.sh <benchmark-dir>"
    exit
fi

E=1e9

BENCHMARK_DIR=$1
NAME=$(basename "$BENCHMARK_DIR")
TARGET=$BENCHMARK_DIR/$NAME.test.fa
QUERY_MSA=$BENCHMARK_DIR/$NAME.train.msa
QUERY_HMM=$BENCHMARK_DIR/$NAME.train.hmm

LONG_SEQ_DIR=$BENCHMARK_DIR/long-seq/
LONG_SEQ_QUERY_DIR=$LONG_SEQ_DIR/query/
LONG_SEQ_TARGET_DIR=$LONG_SEQ_DIR/target/

RESULTS_DIR=$BENCHMARK_DIR/results/nail/

PREP_TIME_1=$RESULTS_DIR/nail.prep.1.time
PREP_TIME_8=$RESULTS_DIR/nail.prep.8.time

SEED_TIME_1=$RESULTS_DIR/nail.seed.1.time
SEED_TIME_8=$RESULTS_DIR/nail.seed.8.time

ALIGN_DEFAULT_1_TIME=$RESULTS_DIR/nail.align.1.default.time
ALIGN_DEFAULT_8_TIME=$RESULTS_DIR/nail.align.8.default.time

ALIGN_FULL_TIME=$RESULTS_DIR/nail.align.full.time

PREP=$RESULTS_DIR/prep/
PREP=$RESULTS_DIR/prep/
SEEDS=$PREP/seeds.json

OUT_DEFAULT=$RESULTS_DIR/nail.default.out
TSV_DEFAULT=$RESULTS_DIR/nail.default.tsv

OUT_FULL=$RESULTS_DIR/nail.full.out
TSV_FULL=$RESULTS_DIR/nail.full.tsv

rm -rf $RESULTS_DIR
mkdir -p $RESULTS_DIR
mkdir $PREP

LONG_SEQ_TSV=$RESULTS_DIR/long-seq.tsv
for ((i=1; i<=6; i++)); do
  LONG_QUERY="$LONG_SEQ_QUERY_DIR${i}.query.fa"
  LONG_TARGET="$LONG_SEQ_TARGET_DIR${i}.target.fa"
  nail search -T tmp.tsv -p $PREP $LONG_QUERY $LONG_TARGET > /dev/null
  cat tmp.tsv >> $LONG_SEQ_TSV
  rm tmp.tsv
done

/usr/bin/time -p -o $PREP_TIME_1 nail prep -t 1 --skip-hmmbuild -p $PREP $QUERY_MSA $TARGET
/usr/bin/time -p -o $PREP_TIME_8 nail prep -t 8 --skip-hmmbuild -p $PREP $QUERY_MSA $TARGET

/usr/bin/time -p -o $SEED_TIME_1 nail seed -t 1 -q $QUERY_HMM -s $SEEDS $PREP
/usr/bin/time -p -o $SEED_TIME_8 nail seed -t 8 -q $QUERY_HMM -s $SEEDS $PREP

/usr/bin/time -p -o $ALIGN_DEFAULT_1_TIME nail align -t 1 -E $E -T $TSV_DEFAULT -O $OUT_DEFAULT $QUERY_HMM $TARGET $SEEDS
/usr/bin/time -p -o $ALIGN_DEFAULT_8_TIME nail align -t 8 -E $E -T $TSV_DEFAULT -O $OUT_DEFAULT $QUERY_HMM $TARGET $SEEDS

/usr/bin/time -p -o $ALIGN_FULL_TIME nail align -t 8 --full-dp -E $E -T $TSV_FULL -O $OUT_FULL $QUERY_HMM $TARGET $SEEDS
