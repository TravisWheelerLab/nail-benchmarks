#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

TMP=./tmp/blast/
mkdir -p $TMP

TARGET_DB=$TMP/target_db
$MAKEBLASTDB -in $TARGET -dbtype prot -out $TARGET_DB > /dev/null

###

run_blastp "blast.seq"
run_psiblast "blast.prf"
