#! /bin/bash

if [ "$#" == 0 ]; then
    echo "usage: ./extract.sh <n>"
    exit
elif [ "$#" == 1 ]; then
    N=$1
fi

shuf pfam.txt > pfam.shuf

LINES="$(head -n $N pfam.shuf)"

while read -r query; do
    # echo "$query"
    esl-afetch pfam.sto $query
done <<< "$LINES"
