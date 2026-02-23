#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

###

set_default E 1e9
set_default RESULTS "${BM_DIR}/results/"

TMP=./tmp/diamond
rm -rf $TMP
mkdir -p $TMP
mkdir -p $RESULTS

TARGET_DB=$TMP/target_db
$DIAMOND makedb --in $TARGET --db $TARGET_DB > /dev/null 2>&1

###

QUERY=$QUERY_FA

# run_diamond "diamond-faster.seq"     "--faster"
# run_diamond "diamond-fast.seq"       "--fast"
# run_diamond "diamond-mid-sens.seq"   "--mid-sensitive"
# run_diamond "diamond-sens.seq"       "--sensitive"
# run_diamond "diamond-more-sens.seq"  "--more-sensitive"
# run_diamond "diamond-very-sens.seq"  "--very-sensitive"

run_diamond "diamond-default.seq"    "--evalue ${E}"
run_diamond "diamond-ultra-sens.seq" "--ultra-sensitive --evalue ${E}"
