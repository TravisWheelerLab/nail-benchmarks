#! /bin/sh

if [ "$#" == 0 ]; then
    echo "usage: ./run-nail.sh <benchmark-dir>"
    exit
fi

E=1e9
BENCHMARK_DIR=$1
FAM_DIR=$BENCHMARK_DIR/fams
TBL=$BENCHMARK_DIR/results.tbl
RESULTS=$BENCHMARK_DIR/nail.tbl

for fam in $FAM_DIR/*; do
  for query in "$fam"/*.q.fa; do
    x=$(basename "$query" .q.fa)
    target=$fam/$x.t.fa
    nail search \
        -t 1 \
        -E $E \
        --tbl-out $TBL \
        --allow-overwrite \
        $query $target
    grep -v ^# $TBL >> $RESULTS
  done
done

