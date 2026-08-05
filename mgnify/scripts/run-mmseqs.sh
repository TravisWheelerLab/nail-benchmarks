#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

set_default THREADS 8
set_default E 10

check_defined S
check_defined BM_DIR
check_defined A
check_defined B
check_defined RESULTS
check_defined TMP

MGY=$BM_DIR/mgy/
mkdir -p $RESULTS

mkdir -p $TMP
ANNOYING=$TMP/annoying

MDB=$TMP/msaDB
$MMSEQS convertmsa $QUERY_MSA $MDB --identifier-field 0 > /dev/null

QDB=$TMP/queryDB-prf
$MMSEQS msa2profile $MDB $QDB --match-mode 1 > /dev/null

TDB=$TMP/targetDB
ADB=$TMP/alignDB

for i in $(seq "$A" "$B"); do
    TARGET="$MGY/$i.fa"
    $MMSEQS createdb $TARGET $TDB > /dev/null
    run_mmseqs "mmseqs.${i}.prf" "-s ${S} --max-seqs 2000 -e ${E}"
done
