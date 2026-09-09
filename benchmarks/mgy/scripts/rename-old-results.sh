#!/usr/bin/env bash
#
# Rename the old shell harness's result tables to the names `mgy import`
# reads, and move them into a pipeline's results directory.
#
# The old harness named a table for the tool and the shard and nothing else --
# `nail.7.prf.tbl` -- and kept the parameters in the directory name
# (`mmseqs-s7.5/`) or hardcoded in the script that ran it (`run-nail.sh` fixes
# `--mmseqs-s 12.0`). So the run name cannot be worked out from the files and
# is given here instead, one directory at a time.
#
# usage: rename-old-results.sh [--dry-run] <src-dir> <run-name> <results-dir>
#
#   src-dir      a directory of <tool>.<shard>.<...>.{tbl,time} files
#   run-name     what the run is called in the new layout, e.g. nail-s12.0
#   results-dir  the pipeline's results/, which is created if missing
#
# Files move rather than copy: a result set worth running elsewhere is not one
# to keep two copies of. Pass --dry-run to see what would happen first.

set -euo pipefail

DRY_RUN=0
if [[ ${1-} == --dry-run ]]; then
    DRY_RUN=1
    shift
fi

if [[ $# -ne 3 ]]; then
    sed -n '/^# usage:/,/^# to keep/p' "$0" >&2
    exit 1
fi

SRC="$1"
RUN="$2"
RESULTS="$3"

if [[ ! -d $SRC ]]; then
    echo "error: no source directory $SRC" >&2
    exit 1
fi

# the run name is a filename component, and a shard is appended to it. one
# holding a path separator would write outside the results directory
if [[ -z $RUN || $RUN == */* ]]; then
    echo "error: run name '$RUN' is a name, not a path" >&2
    exit 1
fi

((DRY_RUN)) || mkdir -p "$RESULTS"

moved=0
skipped=0

for path in "$SRC"/*; do
    [[ -f $path ]] || continue

    file="${path##*/}"

    # <tool>.<shard>.<anything>.<ext>, which is what the old harness wrote.
    # the shard is the second field, the same place the old parser looked for
    # it (mgy-parse.rs, parse_table_indices)
    IFS='.' read -r -a fields <<< "$file"

    if [[ ${#fields[@]} -lt 3 ]]; then
        echo "skip $file (no shard and no extension)" >&2
        ((++skipped))
        continue
    fi

    shard="${fields[1]}"
    ext="${fields[-1]}"

    if [[ ! $shard =~ ^[0-9]+$ ]]; then
        echo "skip $file (second field '$shard' is not a shard number)" >&2
        ((++skipped))
        continue
    fi

    # the seed lists, the stats and the run summaries have no reader in the new
    # layout, and a .time is renamed alongside the table it belongs to
    case $ext in
    tbl | domtbl | time) ;;
    *)
        ((++skipped))
        continue
        ;;
    esac

    dst="$RESULTS/$RUN.$shard.$ext"

    if ((DRY_RUN)); then
        echo "$file -> $RUN.$shard.$ext"
    else
        if [[ -e $dst ]]; then
            echo "error: $dst already exists" >&2
            exit 1
        fi
        mv "$path" "$dst"
    fi

    ((++moved))
done

echo "$moved file(s)$( ((DRY_RUN)) && echo " would move"), $skipped skipped"
