#! /bin/sh

NAME=benchmark

DIR=./output/hmmer/
OUT=$DIR/hmmer.out
TBL=$DIR/hmmer.tbl
DOM=$DIR/hmmer.domtbl

TBL_EVALUE_COL=5

# full sequence E-value
DOMTBL_EVALUE_COL=7

# conditional E-value
# DOMTBL_EVALUE_COL=12

# independent E-value
# DOMTBL_EVALUE_COL=13

mkdir -p $DIR
hmmsearch -E 200 -o $OUT --domtblout $DOM --tblout $TBL ./$NAME/$NAME.hmm ./$NAME/$NAME.fa

grep -v '^#' $TBL | sort -g -k$TBL_EVALUE_COL > $TBL.sorted
