#!/usr/bin/env bash
BLASTP=../tools/bin/blastp
PSIBLAST=../tools/bin/psiblast
MAKEBLASTDB=../tools/bin/makeblastdb

set -e
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

TMP=./tmp/blast/
mkdir -p $TMP
TARGET_DB=$TMP/target_db

TBL_1=$RESULTS/blast.seq.tbl
TIME_1=$RESULTS/blast.seq.time
TBL_2=$RESULTS/blast.cons.tbl
TIME_2=$RESULTS/blast.cons.time
TBL_3=$RESULTS/blast.prf.tbl
TIME_3=$RESULTS/blast.prf.time

$MAKEBLASTDB -in $TARGET -dbtype prot -out $TARGET_DB > /dev/null

echo "running blast seq..."
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
    $BLASTP -query $QUERY_CONS_FA \
    -db $TARGET_DB \
    -out $TBL_2 \
    -outfmt 7 \
    -evalue $E \
    -num_threads $THREADS
cat $TIME_2 | grep real | awk '{print $2 "s"}'
echo

echo "running blast profile..."
start=$(date +%s)
for q in $QUERY_AFA/*.afa; do
  $PSIBLAST -in_msa $q \
  -db $TARGET_DB \
  -outfmt 7 \
  -evalue $E \
  -num_threads $THREADS \
  -comp_based_stats 1 \
  -num_iterations 1 \
  >> $TBL_3
done
end=$(date +%s)
echo "$((end - start))s"
echo "$((end - start))s" >> $TIME_3
