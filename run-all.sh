#! /bin/sh
if [ "$#" == 0 ]; then
    THREADS=8
elif [ "$#" -ge 1 ]; then
    THREADS=$1
fi


./scripts/run-nail.sh $THREADS
./scripts/run-mmseqs.sh $THREADS
#./scripts/run-hmmsearch.sh $THREADS
