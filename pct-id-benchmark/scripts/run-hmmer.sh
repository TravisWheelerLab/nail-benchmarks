#!/usr/bin/env bash
HMMSEARCH=../tools/bin/hmmsearch
PHMMER=../tools/bin/phmmer
HMMBALANCE=../util/scripts/hmmbalance
FASTABALANCE=../util/scripts/fastabalance
TIME_CMD="/usr/bin/time -p"

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

N_SPLITS=$(( THREADS / 4 ))
SPLIT_THREADS=4
SPLIT_DIR=$DIR/query-splits

S_ARGS="--cpu $SPLIT_THREADS -E $E"

TBL=$RESULTS/hmmer.seq.tbl
DOM=$RESULTS/hmmer.seq.domtbl
TIME=$RESULTS/hmmer.seq.time
OUT=$RESULTS/hmmer.seq.out

echo "running phmmer seq..."
SPLIT_TIME=$($TIME_CMD $FASTABALANCE $QUERY_FA $N_SPLITS $SPLIT_DIR 2>&1)
echo "balance time: $(echo $SPLIT_TIME | awk '{print $2 "s"}')"

parallel \
    "${TIME_CMD} -o ${SPLIT_DIR}/{/.}.time \
    $PHMMER $S_ARGS \
    -o ${SPLIT_DIR}/{/.}.out \
    --tblout ${SPLIT_DIR}/{/.}.tbl \
    --domtblout ${SPLIT_DIR}/{/.}.domtbl \
    {} ${TARGET}" \
    ::: "${SPLIT_DIR}"/*.fa

cat $SPLIT_DIR/*.tbl > $TBL
cat $SPLIT_DIR/*.domtbl > $DOM
cat $SPLIT_DIR/*.out > $OUT
cat $SPLIT_DIR/*.time > $TIME

rm -rf $SPLIT_DIR
echo "split times:"
cat $TIME | grep real | awk '{print $2 "s"}'
echo

TBL=$RESULTS/hmmer.prf.tbl
DOM=$RESULTS/hmmer.prf.domtbl
OUT=$RESULTS/hmmer.prf.out
TIME=$RESULTS/hmmer.prf.time

echo "running hmmsearch..."
SPLIT_TIME=$($TIME_CMD $HMMBALANCE $QUERY_HMM $N_SPLITS $SPLIT_DIR 2>&1)
echo "balance time: $(echo $SPLIT_TIME | awk '{print $2 "s"}')"

parallel \
    "${TIME_CMD} -o ${SPLIT_DIR}/{/.}.time \
    $HMMSEARCH $S_ARGS \
    -o ${SPLIT_DIR}/{/.}.out \
    --tblout ${SPLIT_DIR}/{/.}.tbl \
    --domtblout ${SPLIT_DIR}/{/.}.domtbl \
    {} ${TARGET}" \
    ::: "${SPLIT_DIR}"/*.hmm

cat $SPLIT_DIR/*.tbl > $TBL
cat $SPLIT_DIR/*.domtbl > $DOM
cat $SPLIT_DIR/*.out > $OUT
cat $SPLIT_DIR/*.time > $TIME
rm -rf $SPLIT_DIR
echo "split times:"
cat $TIME | grep real | awk '{print $2 "s"}'
echo

# TBL=$RESULTS/hmmer.cons.tbl
# DOM=$RESULTS/hmmer.cons.domtbl
# TIME=$RESULTS/hmmer.cons.time

# echo "running phmmer consensus.."
# SPLIT_TIME=$($TIME_CMD $FASTABALANCE $QUERY_CONS_FA $N_SPLITS $SPLIT_DIR 2>&1)
# echo "balance time: $(echo $SPLIT_TIME | awk '{print $2 "s"}')"

# parallel \
#     "${TIME_CMD} -o ${SPLIT_DIR}/{/.}.time \
#     $PHMMER $S_ARGS \
#     --tblout ${SPLIT_DIR}/{/.}.tbl \
#     --domtblout ${SPLIT_DIR}/{/.}.domtbl \
#     {} ${TARGET}" \
#     ::: "${SPLIT_DIR}"/*.fa

# cat $SPLIT_DIR/*.tbl > $TBL
# cat $SPLIT_DIR/*.domtbl > $DOM
# cat $SPLIT_DIR/*.time > $TIME
# rm -rf $SPLIT_DIR
# echo "split times:"
# cat $TIME | grep real | awk '{print $2 "s"}'
# echo

