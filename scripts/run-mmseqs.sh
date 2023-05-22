#! /bin/sh

NAME=benchmark
TARGET=./$NAME/$NAME.fa
QUERY=./$NAME/$NAME.msa

DIR=./output/mmseqs/
PREP=$DIR/prep/
TARGET_DB=$PREP/targetDb
MSA_DB=$PREP/msaDb
QUERY_DB=$PREP/queryDb
ALIGN_DB=$PREP/alignDb
 
OUT=$DIR/mmseqs.tsv
SORTED=$DIR/results.sorted

EVALUE_COL=7

rm -rf $DIR
mkdir -p $DIR
mkdir $PREP

mmseqs convertmsa $QUERY $MSA_DB --identifier-field 0
mmseqs msa2profile $MSA_DB $QUERY_DB --match-mode 1
mmseqs createdb $TARGET $TARGET_DB
mmseqs search $QUERY_DB $TARGET_DB $ALIGN_DB $PREP -s 7.5 --max-seqs 1000
mmseqs convertalis $QUERY_DB $TARGET_DB $ALIGN_DB $OUT --format-output "target,query,tstart,tend,qstart,qend,evalue"

sort -g -k$EVALUE_COL $OUT > $SORTED