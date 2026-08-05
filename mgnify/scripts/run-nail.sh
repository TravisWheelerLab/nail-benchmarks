#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

set_default THREADS 8
set_default E 10

check_defined BM_DIR
check_defined A
check_defined B
check_defined TMP
check_defined RESULTS

MGY=$BM_DIR/mgy/
mkdir -p $RESULTS

QUERY=$QUERY_HMM
for i in $(seq "$A" "$B"); do
    TARGET="$MGY/$i.fa"
    run_nail "nail.${i}.prf" "--mmseqs-s 12.0 --prog-seed"
done
