#! /bin/sh

if [ "$#" == 0 ]; then
    echo "usage: ./run-mmseqs.sh <benchmark-dir>"
    exit
fi

E=1e9

BENCHMARK_DIR=$1
NAME=$(basename "$BENCHMARK_DIR")
TARGET=$BENCHMARK_DIR/$NAME.test.fa
QUERY=$BENCHMARK_DIR/$NAME.train.msa

RESULTS_DIR=$BENCHMARK_DIR/results/mmseqs/

TIME_DEFAULT_1=$RESULTS_DIR/mmseqs.default.1.time
TIME_DEFAULT_8=$RESULTS_DIR/mmseqs.default.8.time

TIME_SENSITIVE_1=$RESULTS_DIR/mmseqs.sensitive.1.time
TIME_SENSITIVE_8=$RESULTS_DIR/mmseqs.sensitive.8.time

TIME_PREFILTER_NAIL_1=$RESULTS_DIR/mmseqs.prefilter.nail.1.time
TIME_PREFILTER_NAIL_8=$RESULTS_DIR/mmseqs.prefilter.nail.8.time

TIME_ALIGN_NAIL_1=$RESULTS_DIR/mmseqs.align.nail.1.time
TIME_ALIGN_NAIL_8=$RESULTS_DIR/mmseqs.align.nail.8.time

PREP=$RESULTS_DIR/prep/
TARGET_DB=$PREP/targetDb
MSA_DB=$PREP/msaDb
QUERY_DB=$PREP/queryDb

K_NAIL=0
K_SCORE_NAIL=80
MIN_UNGAPPED_SCORE_NAIL=15
MAX_SEQS_NAIL=1000

PREFILTER_DB_NAIL=$PREP/prefilterDb-nail

ALIGN_DB_DEFAULT_1=$PREP/alignDb-default-1
ALIGN_DB_DEFAULT_8=$PREP/alignDb-default-8
ALIGN_DB_SENSITIVE_1=$PREP/alignDb-sensitive-1
ALIGN_DB_SENSITIVE_8=$PREP/alignDb-sensitive-8
ALIGN_DB_NAIL_1=$PREP/alignDb-nail-1
ALIGN_DB_NAIL_8=$PREP/alignDb-nail-8
 
OUT_DEFAULT=$RESULTS_DIR/mmseqs.default.tsv
OUT_SENSITIVE=$RESULTS_DIR/mmseqs.sensitive.tsv
OUT_NAIL=$RESULTS_DIR/mmseqs.nail.tsv

rm -rf $RESULTS_DIR
mkdir -p $RESULTS_DIR
mkdir $PREP

mmseqs convertmsa $QUERY $MSA_DB --identifier-field 0 > /dev/null
mmseqs msa2profile $MSA_DB $QUERY_DB --match-mode 1 > /dev/null
mmseqs createdb $TARGET $TARGET_DB > /dev/null

/usr/bin/time -p -o $TIME_DEFAULT_1 mmseqs search $QUERY_DB $TARGET_DB $ALIGN_DB_DEFAULT_1 $PREP --threads 1 -e $E > /dev/null
/usr/bin/time -p -o $TIME_DEFAULT_8 mmseqs search $QUERY_DB $TARGET_DB $ALIGN_DB_DEFAULT_8 $PREP --threads 8 -e $E > /dev/null
mmseqs convertalis $QUERY_DB $TARGET_DB $ALIGN_DB_DEFAULT_1 $OUT_DEFAULT --format-output "target,query,tstart,tend,qstart,qend,evalue" > /dev/null

/usr/bin/time -p -o $TIME_SENSITIVE_1 mmseqs search $QUERY_DB $TARGET_DB $ALIGN_DB_SENSITIVE_1 $PREP --threads 1 -s 7.5 --max-seqs 1000 -e $E > /dev/null
/usr/bin/time -p -o $TIME_SENSITIVE_8 mmseqs search $QUERY_DB $TARGET_DB $ALIGN_DB_SENSITIVE_8 $PREP --threads 8 -s 7.5 --max-seqs 1000 -e $E > /dev/null
mmseqs convertalis $QUERY_DB $TARGET_DB $ALIGN_DB_SENSITIVE_1 $OUT_SENSITIVE --format-output "target,query,tstart,tend,qstart,qend,evalue" > /dev/null

/usr/bin/time -p -o $TIME_PREFILTER_NAIL_1 mmseqs prefilter $QUERY_DB $TARGET_DB $PREFILTER_DB_NAIL --threads 1 -k $K_NAIL --k-score $K_SCORE_NAIL --min-ungapped-score $MIN_UNGAPPED_SCORE_NAIL --max-seqs $MAX_SEQS_NAIL > /dev/null
/usr/bin/time -p -o $TIME_PREFILTER_NAIL_8 mmseqs prefilter $QUERY_DB $TARGET_DB $PREFILTER_DB_NAIL --threads 8 -k $K_NAIL --k-score $K_SCORE_NAIL --min-ungapped-score $MIN_UNGAPPED_SCORE_NAIL --max-seqs $MAX_SEQS_NAIL > /dev/null
/usr/bin/time -p -o $TIME_ALIGN_NAIL_1 mmseqs align $QUERY_DB $TARGET_DB $PREFILTER_DB_NAIL $ALIGN_DB_NAIL_1 --threads 1 -e $E > /dev/null
/usr/bin/time -p -o $TIME_ALIGN_NAIL_8 mmseqs align $QUERY_DB $TARGET_DB $PREFILTER_DB_NAIL $ALIGN_DB_NAIL_8 --threads 8 -e $E > /dev/null
mmseqs convertalis $QUERY_DB $TARGET_DB $ALIGN_DB_NAIL_1 $OUT_NAIL --format-output "target,query,tstart,tend,qstart,qend,evalue" > /dev/null
