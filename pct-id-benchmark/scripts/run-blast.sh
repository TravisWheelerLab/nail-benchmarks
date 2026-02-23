#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

set_default E 1e9
set_default RESULTS "${BM_DIR}/results/"

TMP=./tmp/blast/
rm -rf $TMP
mkdir -p $TMP
mkdir -p $RESULTS

TARGET_DB=$TMP/target_db
$MAKEBLASTDB -in $TARGET -dbtype prot -out $TARGET_DB > /dev/null

###

run_blastp   "blast.seq" "-evalue ${E}"
run_psiblast "blast.prf" "-evalue ${E}"
