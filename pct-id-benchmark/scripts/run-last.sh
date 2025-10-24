#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

TMP=./tmp/last/
mkdir -p $TMP

TARGET_DB=$TMP/target_db
$LASTDB -p $TARGET_DB $TARGET

###

QUERY=$QUERY_FA
run_last "last.seq"
