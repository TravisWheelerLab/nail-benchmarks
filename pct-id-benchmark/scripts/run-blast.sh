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

$MAKEBLASTDB -in $TARGET -dbtype prot -out $TARGET_DB > /dev/null

TBL_SEQ=$RESULTS/blast.seq.tbl
TIME_SEQ=$RESULTS/blast.seq.time
echo "running blast seq..."
/usr/bin/time -p -o $TIME_SEQ \
    $BLASTP -query $QUERY_FA \
    -db $TARGET_DB \
    -out $TBL_SEQ \
    -outfmt 6 \
    -evalue $E \
    -num_threads $THREADS
cat $TIME_SEQ | grep real | awk '{print $2 "s"}'
echo

TBL_PRF=$RESULTS/blast.prf.tbl
TIME_PRF=$RESULTS/blast.prf.time
echo "running blast profile..."
start=$(date +%s)
for q in $QUERY_AFA/*.afa; do
  $PSIBLAST -in_msa $q \
  -db $TARGET_DB \
  -outfmt 6 \
  -evalue $E \
  -num_threads $THREADS \
  -comp_based_stats 1 \
  -num_iterations 1 \
  >> $TBL_PRF
done
end=$(date +%s)
echo "$((end - start))s"
echo "$((end - start))s" >> $TIME_PRF

# TBL_CONS=$RESULTS/blast.cons.tbl
# TIME_CONS=$RESULTS/blast.cons.time
# echo "running blast consensus..."
# /usr/bin/time -p -o $TIME_CONS \
#     $BLASTP -query $QUERY_CONS_FA \
#     -db $TARGET_DB \
#     -out $TBL_CONS \
#     -outfmt 6 \
#     -evalue $E \
#     -num_threads $THREADS
# cat $TIME_CONS | grep real | awk '{print $2 "s"}'
# echo

