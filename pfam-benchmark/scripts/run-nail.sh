#! /bin/sh

if [ "$#" == 0 ]; then
    echo "usage: ./run-nail.sh <benchmark-dir> [threads]"
    exit
fi

if [ -n "$2" ]; then
    THREADS=$2
else
    THREADS=8
fi

E=1e9

K=6
K_SCORE=80
MIN_UNGAPPED_SCORE=15
MAX_SEQS=2147483647

BENCHMARK_DIR=$1
QUERY=$BENCHMARK_DIR/query.hmm
TARGET=$BENCHMARK_DIR/target.fa

RESULTS_DIR=$BENCHMARK_DIR/results/nail/

TIME_DEFAULT=$RESULTS_DIR/nail.default.time
TIME_DOUBLE=$RESULTS_DIR/nail.double.time
TIME_FULL=$RESULTS_DIR/nail.full.time
TIME_FULL_DOUBLE=$RESULTS_DIR/nail.full-double.time
TIME_NO_FILTERS=$RESULTS_DIR/nail.no-filters.time

TBL_DEFAULT=$RESULTS_DIR/nail.default.tsv
TBL_DOUBLE=$RESULTS_DIR/nail.double.tsv
TBL_FULL=$RESULTS_DIR/nail.full.tsv
TBL_FULL_DOUBLE=$RESULTS_DIR/nail.full-double.tsv
TBL_NO_FILTERS=$RESULTS_DIR/nail.no-filters.tsv

rm -rf $RESULTS_DIR
mkdir -p $RESULTS_DIR
mkdir $PREP

echo "running nail default..."
/usr/bin/time -p -o $TIME_DEFAULT \
    nail search \
    -t $THREADS \
    -E $E \
    --tbl-out $TBL_DEFAULT \
    --mmseqs-k $K \
    --mmseqs-k-score $K_SCORE \
    --mmseqs-min-ungapped-score $MIN_UNGAPPED_SCORE \
    --mmseqs-max-seqs $MAX_SEQS \
    $QUERY $TARGET

awk '/real/ {print "time:", $2}' $TIME_DEFAULT
echo

# echo "running nail double..."
# /usr/bin/time -p -o $TIME_DOUBLE \
#     nail search \
#     -t $THREADS \
#     -E $E \
#     --double-seed \
#     --tbl-out $TBL_DOUBLE \
#     --mmseqs-k $K \
#     --mmseqs-k-score $K_SCORE \
#     --mmseqs-min-ungapped-score $MIN_UNGAPPED_SCORE \
#     --mmseqs-max-seqs $MAX_SEQS \
#     $QUERY $TARGET

# awk '/real/ {print "time:", $2}' $TIME_DOUBLE
# echo

echo "running nail full-dp..."
/usr/bin/time -p -o $TIME_FULL \
    nail search \
    -t $THREADS \
    -E $E \
    --tbl-out $TBL_FULL \
    --full-dp \
    --mmseqs-k $K \
    --mmseqs-k-score $K_SCORE \
    --mmseqs-min-ungapped-score $MIN_UNGAPPED_SCORE \
    --mmseqs-max-seqs $MAX_SEQS \
    $QUERY $TARGET

awk '/real/ {print "time:", $2}' $TIME_FULL
echo

# echo "running nail full-dp double..."
# /usr/bin/time -p -o $TIME_FULL_DOUBLE \
#     nail search \
#     -t $THREADS \
#     -E $E \
#     --tbl-out $TBL_FULL_DOUBLE \
#     --double-seed \
#     --full-dp \
#     --mmseqs-k $K \
#     --mmseqs-k-score $K_SCORE \
#     --mmseqs-min-ungapped-score $MIN_UNGAPPED_SCORE \
#     --mmseqs-max-seqs $MAX_SEQS \
#     $QUERY $TARGET

# awk '/real/ {print "time:", $2}' $TIME_FULL_DOUBLE
# echo

echo "running nail no-filters..."
/usr/bin/time -p -o $TIME_NO_FILTERS \
    nail search \
    -t $THREADS \
    -E $E \
    --tbl-out $TBL_NO_FILTERS \
    -F 1e9 \
    -C 1e9 \
    --mmseqs-k $K \
    --mmseqs-k-score $K_SCORE \
    --mmseqs-min-ungapped-score $MIN_UNGAPPED_SCORE \
    --mmseqs-max-seqs $MAX_SEQS \
    $QUERY $TARGET

awk '/real/ {print "time:", $2}' $TIME_NO_FILTERS
echo
