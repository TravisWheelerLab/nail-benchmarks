#! /bin/bash

if [ "$#" == 0 ]; then
    echo "usage: ./extract.sh <n>"
    exit
elif [ "$#" == 1 ]; then
    N=$1
fi

shuf pfam.txt | head -n $N > tmp

# while read -r query; do
#     esl-afetch pfam.sto $query
# done <<< "$LINES"

esl-afetch -f pfam.sto tmp
