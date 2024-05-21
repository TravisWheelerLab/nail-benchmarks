#! /bin/sh

if [ "$#" == 0 ]; then
    echo "usage: ./run-mmseqs.sh <benchmark-dir> [threads]"
    exit
fi


if [ -n "$2" ]; then
    THREADS=$2
else
    THREADS=1
fi

E=1e9

BENCHMARK_DIR=$1
NAME=$(basename "$BENCHMARK_DIR")
TARGET=$BENCHMARK_DIR/$NAME.test.fa
QUERY=$BENCHMARK_DIR/$NAME.train.msa

RESULTS_DIR=$BENCHMARK_DIR/results/mmseqs/

TIME_DEFAULT=$RESULTS_DIR/mmseqs.default.time
TIME_SENSITIVE=$RESULTS_DIR/mmseqs.sensitive.time
TIME_NAIL=$RESULTS_DIR/mmseqs.nail.time

PREP=$RESULTS_DIR/prep/
TARGET_DB=$PREP/targetDb
MSA_DB=$PREP/msaDb
QUERY_DB=$PREP/queryDb

# K_DEFAULT=6
# K_SENSITIVE=6
# MAX_SEQS_SENSITIVE=1000

K_NAIL=6
K_SCORE_NAIL=80
MIN_UNGAPPED_SCORE_NAIL=15
MAX_SEQS_NAIL=1000

ALIGN_DB_DEFAULT=$PREP/alignDb-default
ALIGN_DB_SENSITIVE=$PREP/alignDb-sensitive
ALIGN_DB_NAIL=$PREP/alignDb-nail
 
OUT_DEFAULT=$RESULTS_DIR/mmseqs.default.tsv
OUT_SENSITIVE=$RESULTS_DIR/mmseqs.sensitive.tsv
OUT_NAIL=$RESULTS_DIR/mmseqs.nail.tsv

STDOUT=$RESULTS_DIR/tmp

rm -rf $RESULTS_DIR
mkdir -p $RESULTS_DIR
mkdir $PREP

echo "preparing files..."
mmseqs convertmsa $QUERY $MSA_DB \
    --identifier-field 0 \
    > /dev/null

mmseqs msa2profile $MSA_DB $QUERY_DB \
    --match-mode 1 \
    > /dev/null

mmseqs createdb $TARGET $TARGET_DB \
    > /dev/null

echo "running mmseqs default..."
/usr/bin/time -p -o $TIME_DEFAULT \
    mmseqs search $QUERY_DB $TARGET_DB $ALIGN_DB_DEFAULT $PREP \
    --threads $THREADS \
    -e $E \
    > $STDOUT

awk '/real/ {print "time:", $2}' $TIME_DEFAULT
grep "Index table k-mer threshold:" $STDOUT
grep "k-mer similarity threshold: " $STDOUT
echo

mmseqs convertalis $QUERY_DB $TARGET_DB $ALIGN_DB_DEFAULT $OUT_DEFAULT \
    --format-output "target,query,tstart,tend,qstart,qend,evalue" \
    > /dev/null

echo "running mmseqs sensitive..."
/usr/bin/time -p -o $TIME_SENSITIVE \
    mmseqs search $QUERY_DB $TARGET_DB $ALIGN_DB_SENSITIVE $PREP \
    --threads $THREADS \
    -s 7.5 \
    --max-seqs 1000 \
    -e $E \
    > $STDOUT

awk '/real/ {print "time:", $2}' $TIME_SENSITIVE
grep "Index table k-mer threshold:" $STDOUT
grep "k-mer similarity threshold: " $STDOUT
echo

mmseqs convertalis $QUERY_DB $TARGET_DB $ALIGN_DB_SENSITIVE $OUT_SENSITIVE \
    --format-output "target,query,tstart,tend,qstart,qend,evalue" \
    > /dev/null

echo "running mmseqs nail pipeline settings..."
/usr/bin/time -p -o $TIME_NAIL \
    mmseqs search $QUERY_DB $TARGET_DB $ALIGN_DB_NAIL $PREP \
    --threads $THREADS \
    -k $K_NAIL \
    --k-score $K_SCORE_NAIL \
    --min-ungapped-score $MIN_UNGAPPED_SCORE_NAIL \
    --max-seqs $MAX_SEQS_NAIL \
    -e $E \
    > $STDOUT

awk '/real/ {print "time:", $2}' $TIME_NAIL
grep "Index table k-mer threshold:" $STDOUT
grep "k-mer similarity threshold: " $STDOUT
echo

mmseqs convertalis $QUERY_DB $TARGET_DB $ALIGN_DB_NAIL $OUT_NAIL \
    --format-output "target,query,tstart,tend,qstart,qend,evalue" \
    > /dev/null

sort -k7g -o $OUT_DEFAULT $OUT_DEFAULT
sort -k7g -o $OUT_SENSITIVE $OUT_SENSITIVE
sort -k7g -o $OUT_NAIL $OUT_NAIL
