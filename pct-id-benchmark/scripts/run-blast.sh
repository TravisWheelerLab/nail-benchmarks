#!/usr/bin/env bash
BLASTP=../tools/bin/blastp
MAKEBLASTDB=../tools/bin/makeblastdb

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

TMP=./tmp/blast/
mkdir -p $TMP
TARGET_DB=$TMP/target_db

TBL_1=$RESULTS/blast.pair.tbl
TIME_1=$RESULTS/blast.pair.time
TBL_2=$RESULTS/blast.cons.tbl
TIME_2=$RESULTS/blast.cons.time

$MAKEBLASTDB -in $TARGET -dbtype prot -out $TARGET_DB > /dev/null

echo "running blast pairwise..."
/usr/bin/time -p -o $TIME_1 \
    $BLASTP -query $QUERY_FA \
    -db $TARGET_DB \
    -out $TBL_1 \
    -outfmt 7 \
    -evalue $E \
    -num_threads $THREADS
cat $TIME_1 | grep real | awk '{print $2 "s"}'
echo

echo "running blast consensus..."
/usr/bin/time -p -o $TIME_2 \
    $BLASTP -query $QUERY_FA \
    -db $TARGET_DB \
    -out $TBL_2 \
    -outfmt 7 \
    -evalue $E \
    -num_threads $THREADS
cat $TIME_2 | grep real | awk '{print $2 "s"}'
echo
