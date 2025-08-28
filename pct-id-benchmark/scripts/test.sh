#!/usr/bin/env bash
DIR="$(cd -- "$(dirname -- "$0")" && pwd)"

. $DIR/vars.sh
set_vars "$@"

if [ "$THREADS" -lt 10 ] || [ "$THREADS" -gt 50 ]; then
  echo "bin: $THREADS"
  echo "usage: ./test.sh <benchmark-dir> [10<=bin<=50]"
  exit
else
    BIN=$THREADS
fi

if [ ! -f "$TARGET.ssi" ]; then
    esl-sfetch --index $TARGET
fi

TMP=./tmp/
mkdir -p $TMP

TMP_Q=$TMP/q.hmm
TMP_T=$TMP/t.fa

TMP_DOM=$TMP/tmp.domtbl
TMP_TBL=$TMP/tmp.tbl

grep $BIN% $BENCH_TBL | while read -r target domain query; do
  [[ $target =~ ^# ]] && continue

  fam=$(echo "$target" | awk -F'|' '{print $1}')
  range=$(echo "$target" | awk -F'|' '{print $2}')
  s=$(echo "$range" | awk -F'-' '{print $1}')
  e=$(echo "$range" | awk -F'-' '{print $2}')
  l=$((e - s))

  echo $l
  # echo ">($target) (T:$domain) (Q:$query)"
  # hmmfetch $QUERY_HMM $fam > $TMP_Q
  # esl-sfetch $TARGET $target > $TMP_T

  # hmmsearch \
  #   --domtblout $TMP_DOM \
  #   --tblout $TMP_TBL \
  #   $TMP_Q $TMP_T > /dev/null

  # res=$(grep -v '^#' $TMP_TBL | awk '{print $5, $6, $8, $9}')
  # if [ -z "$res" ]; then
  #   echo "[NONE]"
  # else
  #     echo "$res"
  # fi

  # echo


done
