#! /bin/sh

if [ "$#" == 0 ]; then
    THREADS=8
elif [ "$#" -ge 1 ]; then
    THREADS=$1
fi

NAME=benchmark
TARGET=./$NAME/$NAME.fa
QUERY=./$NAME/$NAME.msa

DIR=./output/mmseqs/
TIME=$DIR/mmseqs.time

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

mmseqs convertmsa $QUERY $MSA_DB --identifier-field 0 > /dev/null
mmseqs msa2profile $MSA_DB $QUERY_DB --match-mode 1 > /dev/null
mmseqs createdb $TARGET $TARGET_DB > /dev/null
/usr/bin/time -ph -o $TIME mmseqs search $QUERY_DB $TARGET_DB $ALIGN_DB $PREP --threads $THREADS -s 7.5 --max-seqs 1000 > /dev/null
mmseqs convertalis $QUERY_DB $TARGET_DB $ALIGN_DB $OUT --format-output "target,query,tstart,tend,qstart,qend,evalue" > /dev/null

sort -g -k$EVALUE_COL $OUT > $SORTED