#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

set_default E 1e9
set_default RESULTS "${BM_DIR}/results/"

TMP=./tmp/last/
rm -rf $TMP
mkdir -p $TMP
mkdir -p $RESULTS

TARGET_DB=$TMP/target_db
$LASTDB -p $TARGET_DB $TARGET

###

QUERY=$QUERY_FA
run_last "last.seq" "-E ${E}"
