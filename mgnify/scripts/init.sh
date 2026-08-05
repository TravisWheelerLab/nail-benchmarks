#!/usr/bin/env bash

init() {
    set -e

    local SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
    local UTIL_DIR="$(realpath $SCRIPT_DIR/../../util)"
    . "$UTIL_DIR/scripts/run.sh"

    parse_args "$@"
    set_numa_prefix
    set_tool_vars

    if [ "${#POS_ARGS[@]}" -ne 1 ]; then
        echo "usage: $0 [options] <benchmark_dir>" >&2
        exit 1
    fi
    
    BM_DIR="${POS_ARGS[0]}"
    if [ -d "$BM_DIR" ]; then
        echo "benchmark: $BM_DIR"
        BM_DIR="$(realpath "${POS_ARGS[0]}")"
    else
        echo "error: benchmark directory $BM_DIR invalid"
        exit 1
    fi

    QUERY_HMM="$BM_DIR/query.hmm"
    QUERY_MSA="$BM_DIR/query.sto"
    TARGET="$BM_DIR/target.fa"
}

init "$@"
