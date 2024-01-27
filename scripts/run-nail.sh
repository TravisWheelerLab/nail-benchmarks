#! /bin/sh

if [ "$#" == 0 ]; then
    echo "usage: ./run-nail.sh <benchmark-dir>"
    exit
fi

E=1e9
THREADS=8

K=6
K_SCORE=80
MIN_UNGAPPED_SCORE=15
MAX_SEQS=1000

BENCHMARK_DIR=$1
NAME=$(basename "$BENCHMARK_DIR")
TARGET=$BENCHMARK_DIR/$NAME.test.fa
QUERY_MSA=$BENCHMARK_DIR/$NAME.train.msa
QUERY_HMM=$BENCHMARK_DIR/$NAME.train.hmm

LONG_SEQ_DIR=$BENCHMARK_DIR/long-seq/
LONG_SEQ_QUERY_DIR=$LONG_SEQ_DIR/query/
LONG_SEQ_TARGET_DIR=$LONG_SEQ_DIR/target/

RESULTS_DIR=$BENCHMARK_DIR/results/nail/

PREP_TIME=$RESULTS_DIR/nail.prep.time
SEED_TIME=$RESULTS_DIR/nail.seed.time
ALIGN_DEFAULT_TIME=$RESULTS_DIR/nail.align.default.time
ALIGN_FULL_TIME=$RESULTS_DIR/nail.align.full.time
ALIGN_NO_FILTERS_TIME=$RESULTS_DIR/nail.align.no-filters.time

PREP=$RESULTS_DIR/prep/
PREP=$RESULTS_DIR/prep/
SEEDS=$PREP/seeds.json

OUT_DEFAULT=$RESULTS_DIR/nail.default.out
TSV_DEFAULT=$RESULTS_DIR/nail.default.tsv

OUT_FULL=$RESULTS_DIR/nail.full.out
TSV_FULL=$RESULTS_DIR/nail.full.tsv

OUT_FULL=$RESULTS_DIR/nail.full.out
TSV_FULL=$RESULTS_DIR/nail.full.tsv

OUT_NO_FILTERS=$RESULTS_DIR/nail.no-filters.out
TSV_NO_FILTERS=$RESULTS_DIR/nail.no-filters.tsv

rm -rf $RESULTS_DIR
mkdir -p $RESULTS_DIR
mkdir $PREP

echo "running nail on long sequence pairs..."
LONG_SEQ_TSV=$RESULTS_DIR/long-seq.tsv
for ((i=1; i<=6; i++)); do
  LONG_QUERY="$LONG_SEQ_QUERY_DIR${i}.query.fa"
  LONG_TARGET="$LONG_SEQ_TARGET_DIR${i}.target.fa"
  nail search -T tmp.tsv -p $PREP $LONG_QUERY $LONG_TARGET > /dev/null
  cat tmp.tsv >> $LONG_SEQ_TSV
  rm tmp.tsv
done

echo "running nail prep..."
/usr/bin/time -p -o $PREP_TIME \
    nail prep \
    -t $THREADS \
    --skip-hmmbuild \
    -p $PREP \
    $QUERY_MSA $TARGET

awk '/real/ {print "time:", $2}' $PREP_TIME
echo

echo "running nail seed..."
/usr/bin/time -p -o $SEED_TIME \
    nail seed \
    -t $THREADS \
    -q $QUERY_HMM \
    -s $SEEDS \
    --mmseqs_k $K \
    --mmseqs_k_score $K_SCORE \
    --mmseqs_min_ungapped_score $MIN_UNGAPPED_SCORE \
    --mmseqs_max_seqs $MAX_SEQS \
    $PREP

awk '/real/ {print "time:", $2}' $SEED_TIME
echo

echo "running nail align default..."
/usr/bin/time -p -o $ALIGN_DEFAULT_TIME \
    nail align \
    -t $THREADS \
    -E $E \
    -T $TSV_DEFAULT \
    -O $OUT_DEFAULT \
    $QUERY_HMM $TARGET $SEEDS

awk '/real/ {print "time:", $2}' $ALIGN_DEFAULT_TIME
echo

echo "running nail align full-dp..."
/usr/bin/time -p -o $ALIGN_FULL_TIME \
    nail align \
    -t $THREADS \
    --full-dp \
    -E $E \
    -T $TSV_FULL \
    -O $OUT_FULL \
    $QUERY_HMM $TARGET $SEEDS

awk '/real/ {print "time:", $2}' $ALIGN_FULL_TIME
echo

echo "running nail align no-filters..."
/usr/bin/time -p -o $ALIGN_NO_FILTERS_TIME \
    nail align \
    -t $THREADS \
    --forward-thresh 1e9 \
    --cloud-thresh 1e9 \
    -E $E \
    -T $TSV_NO_FILTERS \
    -O $OUT_NO_FILTERS \
    $QUERY_HMM $TARGET $SEEDS

awk '/real/ {print "time:", $2}' $ALIGN_NO_FILTERS_TIME
echo
