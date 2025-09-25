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

S_ARGS="--cpu $SPLIT_THREADS -E $E -o /dev/null"

TBL_SEQ=$RESULTS/hmmer.seq.tbl
DOM_SEQ=$RESULTS/hmmer.seq.domtbl
TIME_SEQ=$RESULTS/hmmer.seq.time

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

cat $SPLIT_DIR/*.tbl > $TBL_SEQ
cat $SPLIT_DIR/*.domtbl > $DOM_SEQ
cat $SPLIT_DIR/*.time > $TIME_SEQ
rm -rf $SPLIT_DIR
echo "split times:"
cat $TIME_SEQ | grep real | awk '{print $2 "s"}'
echo

TBL_PRF=$RESULTS/hmmer.prf.tbl
DOM_PRF=$RESULTS/hmmer.prf.domtbl
TIME_PRF=$RESULTS/hmmer.prf.time

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

cat $SPLIT_DIR/*.tbl > $TBL_PRF
cat $SPLIT_DIR/*.domtbl > $DOM_PRF
cat $SPLIT_DIR/*.time > $TIME_PRF
rm -rf $SPLIT_DIR
echo "split times:"
cat $TIME_PRF | grep real | awk '{print $2 "s"}'
echo

# TBL_CONS=$RESULTS/hmmer.cons.tbl
# DOM_CONS=$RESULTS/hmmer.cons.domtbl
# TIME_CONS=$RESULTS/hmmer.cons.time

# echo "running phmmer consensus.."
# SPLIT_TIME=$($TIME $FASTABALANCE $QUERY_CONS_FA $N_SPLITS $SPLIT_DIR 2>&1)
# echo "balance time: $(echo $SPLIT_TIME | awk '{print $2 "s"}')"

# parallel \
#     "${TIME} -o ${SPLIT_DIR}/{/.}.time \
#     $PHMMER $S_ARGS \
#     --tblout ${SPLIT_DIR}/{/.}.tbl \
#     --domtblout ${SPLIT_DIR}/{/.}.domtbl \
#     {} ${TARGET}" \
#     ::: "${SPLIT_DIR}"/*.fa

# cat $SPLIT_DIR/*.tbl > $TBL_CONS
# cat $SPLIT_DIR/*.domtbl > $DOM_CONS
# cat $SPLIT_DIR/*.time > $TIME_CONS
# rm -rf $SPLIT_DIR
# echo "split times:"
# cat $TIME_CONS | grep real | awk '{print $2 "s"}'
# echo

