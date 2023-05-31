#! /bin/sh

if [ "$#" == 0 ]; then
    THREADS=8
elif [ "$#" -ge 1 ]; then
    THREADS=$1
fi

NAME=benchmark
TARGET=./$NAME/$NAME.fa
QUERY=./$NAME/$NAME.hmm

DIR=./output/hmmer/
TIME=$DIR/hmmer.time
OUT=$DIR/hmmer.out
TBL=$DIR/hmmer.tbl
DOM=$DIR/hmmer.domtbl
SORTED=$DIR/results.sorted

TBL_EVALUE_COL=5

# full sequence E-value
DOMTBL_EVALUE_COL=7

# conditional E-value
# DOMTBL_EVALUE_COL=12

# independent E-value
# DOMTBL_EVALUE_COL=13

mkdir -p $DIR
/usr/bin/time -p -o $TIME hmmsearch --cpu $THREADS -E 200 -o $OUT --domtblout $DOM --tblout $TBL $QUERY $TARGET

grep -v '^#' $TBL | sort -g -k$TBL_EVALUE_COL > $SORTED


