#! /bin/sh
HMMSEARCH=../tools/bin/hmmsearch
PHMMER=../tools/bin/phmmer

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

O1=$RESULTS/hmmer.seq.tbl
T1=$RESULTS/hmmer.seq.time
O2=$RESULTS/hmmer.hmm.tbl
T2=$RESULTS/hmmer.hmm.time

S_ARGS="--cpu $THREADS -E $E"

echo "running phmmer.."
/usr/bin/time -p -o $T1 \
    $PHMMER $S_ARGS \
    --tblout $O1 \
    $QUERY_FA $TARGET > /dev/null
cat $T1 | grep real

echo "running hmmsearch ..."
/usr/bin/time -p -o $T2 \
    $HMMSEARCH $S_ARGS \
    --tblout $O2 \
     $QUERY_HMM $TARGET > /dev/null
cat $T2 | grep real
