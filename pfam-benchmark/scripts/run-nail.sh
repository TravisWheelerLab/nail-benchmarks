#! /bin/sh

if [ "$#" == 0 ]; then
    echo "usage: ./run-nail.sh <benchmark-dir> [threads]"
    exit
fi

if [ -n "$2" ]; then
    THREADS=$2
else
    THREADS=1
fi

E=1e9

K=6
K_SCORE=80
MIN_UNGAPPED_SCORE=15
MAX_SEQS=1000

BENCHMARK_DIR=$1
NAME=$(basename "$BENCHMARK_DIR")
TARGET=$BENCHMARK_DIR/$NAME.test.fa
# QUERY=$BENCHMARK_DIR/$NAME.train.msa
QUERY=$BENCHMARK_DIR/$NAME.train.hmm

LONG_SEQ_DIR=$BENCHMARK_DIR/long-seq/
LONG_SEQ_QUERY_DIR=$LONG_SEQ_DIR/query/
LONG_SEQ_TARGET_DIR=$LONG_SEQ_DIR/target/

RESULTS_DIR=$BENCHMARK_DIR/results/nail/

PREP=$RESULTS_DIR/prep/

TIME_DEFAULT=$RESULTS_DIR/nail.default.time
OUT_DEFAULT=$RESULTS_DIR/nail.default.out
TSV_DEFAULT=$RESULTS_DIR/nail.default.tsv

TIME_FULL=$RESULTS_DIR/nail.full.time
OUT_FULL=$RESULTS_DIR/nail.full.out
TSV_FULL=$RESULTS_DIR/nail.full.tsv

TIME_NO_FILTERS=$RESULTS_DIR/nail.no-filters.time
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
  nail search -T tmp.tsv --prep $PREP $LONG_QUERY $LONG_TARGET > /dev/null
  cat tmp.tsv >> $LONG_SEQ_TSV
  rm tmp.tsv
done

echo "running nail default..."
/usr/bin/time -p -o $TIME_DEFAULT \
    nail search \
    -t $THREADS \
    -E $E \
    -T $TSV_DEFAULT \
    -O $OUT_DEFAULT \
    --prep $PREP \
    --mmseqs-k $K \
    --mmseqs-k-score $K_SCORE \
    --mmseqs-min-ungapped_score $MIN_UNGAPPED_SCORE \
    --mmseqs-max-seqs $MAX_SEQS \
    $QUERY $TARGET

awk '/real/ {print "time:", $2}' $TIME_DEFAULT
echo

echo "running nail full-dp..."
/usr/bin/time -p -o $TIME_FULL \
    nail search \
    -t $THREADS \
    -E $E \
    -T $TSV_FULL \
    -O $OUT_FULL \
    --prep $PREP \
    --full-dp \
    --mmseqs-k $K \
    --mmseqs-k-score $K_SCORE \
    --mmseqs-min-ungapped_score $MIN_UNGAPPED_SCORE \
    --mmseqs-max-seqs $MAX_SEQS \
    $QUERY $TARGET

awk '/real/ {print "time:", $2}' $TIME_FULL
echo

echo "running nail no-filters..."
/usr/bin/time -p -o $TIME_NO_FILTERS \
    nail search \
    -t $THREADS \
    -E $E \
    -T $TSV_NO_FILTERS \
    -O $OUT_NO_FILTERS \
    --prep $PREP \
    --forward-thresh 1e9 \
    --cloud-thresh 1e9 \
    --mmseqs-k $K \
    --mmseqs-k-score $K_SCORE \
    --mmseqs-min-ungapped_score $MIN_UNGAPPED_SCORE \
    --mmseqs-max-seqs $MAX_SEQS \
    $QUERY $TARGET

awk '/real/ {print "time:", $2}' $TIME_NO_FILTERS
echo
