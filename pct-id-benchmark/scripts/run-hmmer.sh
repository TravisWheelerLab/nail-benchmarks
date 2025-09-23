#!/usr/bin/env bash
HMMSEARCH=../tools/bin/hmmsearch
PHMMER=../tools/bin/phmmer
HMMBALANCE=../util/scripts/hmmbalance
FASTABALANCE=../util/scripts/fastabalance
TIME="/usr/bin/time -p"

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

N_SPLITS=$(( THREADS / 4 ))
SPLIT_THREADS=4
SPLIT_DIR=$DIR/query-splits

TBL_1=$RESULTS/hmmer.seq.tbl
DOM_1=$RESULTS/hmmer.seq.domtbl
TIME_1=$RESULTS/hmmer.seq.time

TBL_2=$RESULTS/hmmer.cons.tbl
DOM_2=$RESULTS/hmmer.cons.domtbl
TIME_2=$RESULTS/hmmer.cons.time

TBL_3=$RESULTS/hmmer.prf.tbl
DOM_3=$RESULTS/hmmer.prf.domtbl
TIME_3=$RESULTS/hmmer.prf.time

S_ARGS="--cpu $SPLIT_THREADS -E $E -o /dev/null"

echo "running phmmer seq..."
SPLIT_TIME=$($TIME $FASTABALANCE $QUERY_FA $N_SPLITS $SPLIT_DIR 2>&1)
echo "balance time: $(echo $SPLIT_TIME | awk '{print $2 "s"}')"

parallel \
    "${TIME} -o ${SPLIT_DIR}/{/.}.time \
    $PHMMER $S_ARGS \
    --tblout ${SPLIT_DIR}/{/.}.tbl \
    --domtblout ${SPLIT_DIR}/{/.}.domtbl \
    {} ${TARGET}" \
    ::: "${SPLIT_DIR}"/*.fa

cat $SPLIT_DIR/*.tbl > $TBL_1
cat $SPLIT_DIR/*.domtbl > $DOM_1
cat $SPLIT_DIR/*.time > $TIME_1
rm -rf $SPLIT_DIR
echo "split times:"
cat $TIME_1 | grep real | awk '{print $2 "s"}'
echo

echo "running phmmer consensus.."
SPLIT_TIME=$($TIME $FASTABALANCE $QUERY_CONS_FA $N_SPLITS $SPLIT_DIR 2>&1)
echo "balance time: $(echo $SPLIT_TIME | awk '{print $2 "s"}')"

parallel \
    "${TIME} -o ${SPLIT_DIR}/{/.}.time \
    $PHMMER $S_ARGS \
    --tblout ${SPLIT_DIR}/{/.}.tbl \
    --domtblout ${SPLIT_DIR}/{/.}.domtbl \
    {} ${TARGET}" \
    ::: "${SPLIT_DIR}"/*.fa

cat $SPLIT_DIR/*.tbl > $TBL_2
cat $SPLIT_DIR/*.domtbl > $DOM_2
cat $SPLIT_DIR/*.time > $TIME_2
rm -rf $SPLIT_DIR
echo "split times:"
cat $TIME_2 | grep real | awk '{print $2 "s"}'
echo

echo "running hmmsearch..."
SPLIT_TIME=$($TIME $HMMBALANCE $QUERY_HMM $N_SPLITS $SPLIT_DIR 2>&1)
echo "balance time: $(echo $SPLIT_TIME | awk '{print $2 "s"}')"

parallel \
    "${TIME} -o ${SPLIT_DIR}/{/.}.time \
    $HMMSEARCH $S_ARGS \
    --tblout ${SPLIT_DIR}/{/.}.tbl \
    --domtblout ${SPLIT_DIR}/{/.}.domtbl \
    {} ${TARGET}" \
    ::: "${SPLIT_DIR}"/*.hmm

cat $SPLIT_DIR/*.tbl > $TBL_3
cat $SPLIT_DIR/*.domtbl > $DOM_3
cat $SPLIT_DIR/*.time > $TIME_3
rm -rf $SPLIT_DIR
echo "split times:"
cat $TIME_3 | grep real | awk '{print $2 "s"}'
echo

