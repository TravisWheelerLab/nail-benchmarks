#! /bin/sh

if [ "$#" == 0 ]; then
    echo "usage: ./run-mmseqs.sh <benchmark-dir> [threads]"
    exit
elif [ "$#" == 1 ]; then
    THREADS=8
elif [ "$#" -ge 2 ]; then
    THREADS=$2
fi

MMSEQS=mmseqs

BENCHMARK_DIR=$1
NAME=$(basename "$BENCHMARK_DIR")
TARGET=$BENCHMARK_DIR/$NAME.test.fa
QUERY=$BENCHMARK_DIR/$NAME.train.msa

RESULTS_DIR=$BENCHMARK_DIR/results/mmseqs/

TIME_DEFAULT=$RESULTS_DIR/mmseqs.default.time
TIME_SENSITIVE=$RESULTS_DIR/mmseqs.sensitive.time
TIME_PREFILTER_NAIL=$RESULTS_DIR/mmseqs.prefilter.nail.time
TIME_ALIGN_NAIL=$RESULTS_DIR/mmseqs.align.nail.time

PREP=$RESULTS_DIR/prep/
TARGET_DB=$PREP/targetDb
MSA_DB=$PREP/msaDb
QUERY_DB=$PREP/queryDb

K_NAIL=0
K_SCORE_NAIL=80
MIN_UNGAPPED_SCORE_NAIL=15
MAX_SEQS_NAIL=1000

PREFILTER_DB_NAIL=$PREP/prefilterDb-nail

ALIGN_DB_DEFAULT=$PREP/alignDb-default
ALIGN_DB_SENSITIVE=$PREP/alignDb-sensitive
ALIGN_DB_NAIL=$PREP/alignDb-nail
 
OUT_DEFAULT=$RESULTS_DIR/mmseqs.default.tsv
OUT_SENSITIVE=$RESULTS_DIR/mmseqs.sensitive.tsv
OUT_NAIL=$RESULTS_DIR/mmseqs.nail.tsv

rm -rf $RESULTS_DIR
mkdir -p $RESULTS_DIR
mkdir $PREP

$MMSEQS convertmsa $QUERY $MSA_DB --identifier-field 0 > /dev/null
$MMSEQS msa2profile $MSA_DB $QUERY_DB --match-mode 1 > /dev/null
$MMSEQS createdb $TARGET $TARGET_DB > /dev/null

/usr/bin/time -p -o $TIME_DEFAULT $MMSEQS search $QUERY_DB $TARGET_DB $ALIGN_DB_DEFAULT $PREP --threads $THREADS -e 200 > /dev/null
$MMSEQS convertalis $QUERY_DB $TARGET_DB $ALIGN_DB_DEFAULT $OUT_DEFAULT --format-output "target,query,tstart,tend,qstart,qend,evalue" > /dev/null

/usr/bin/time -p -o $TIME_SENSITIVE $MMSEQS search $QUERY_DB $TARGET_DB $ALIGN_DB_SENSITIVE $PREP --threads $THREADS -s 7.5 --max-seqs 1000 -e 200 > /dev/null
$MMSEQS convertalis $QUERY_DB $TARGET_DB $ALIGN_DB_SENSITIVE $OUT_SENSITIVE --format-output "target,query,tstart,tend,qstart,qend,evalue" > /dev/null

/usr/bin/time -p -o $TIME_PREFILTER_NAIL $MMSEQS prefilter $QUERY_DB $TARGET_DB $PREFILTER_DB_NAIL --threads $THREADS -k $K_NAIL --k-score $K_SCORE_NAIL --min-ungapped-score $MIN_UNGAPPED_SCORE_NAIL --max-seqs $MAX_SEQS_NAIL > /dev/null
/usr/bin/time -p -o $TIME_ALIGN_NAIL $MMSEQS align $QUERY_DB $TARGET_DB $PREFILTER_DB_NAIL $ALIGN_DB_NAIL --threads $THREADS -e 200 > /dev/null
$MMSEQS convertalis $QUERY_DB $TARGET_DB $ALIGN_DB_NAIL $OUT_NAIL --format-output "target,query,tstart,tend,qstart,qend,evalue" > /dev/null
