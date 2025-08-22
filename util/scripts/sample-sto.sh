#! /bin/bash

if [ "$#" == 0 ]; then
    echo "usage: ./sample-sto.sh <x.sto> <n>"
    exit
elif [ "$#" == 2 ]; then
    F=$1
    N=$2
fi

esl-alistat $F | grep name | awk '{print $3}' > tmp.names

shuf tmp.names | head -n $N > tmp.names.shuf

esl-afetch -f $F tmp.names.shuf

rm tmp.names 
rm tmp.names.shuf
