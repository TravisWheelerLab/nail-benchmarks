check_defined() {
    # the arg here needs to be the string name 
    # of the variable, not the variable itself
    local VAR=$1
    local CTX=$2 
    
    # in bash, ${!VAR} means expand the 
    # variable whose name is stored in VAR
    if [[ -z "${!VAR}" ]]; then
        echo "error: var $VAR undefined" >&2
        [[ -n "$CTX" ]] && echo "context: $CTX"

        exit 1
    fi
}

set_default() {
    local VAR="$1" 
    local VAL="$2"
    if [[ -z "${!VAR}" ]]; then
        printf -v "$VAR" '%s' "$VAL"
    fi
}

parse_args() {
    POS_ARGS=()
    NAMED_ARGS=()
    while [[ $# -gt 0 ]]; do
        if [[ $1 == --* && $# -gt 1 ]]; then
            local var="${1#--}"
            var=$(printf '%s' "$var" | tr '[:lower:]' '[:upper:]')
            local val="$2"
            printf -v "$var" '%s' "$val"
            NAMED_ARGS+=("$var")
            shift 2
        else
            POS_ARGS+=("$1")
            shift
        fi
    done
}

print_args() {
    for v in "${NAMED_ARGS[@]}"; do
        echo "$v=${!v}"
    done

    echo "positional:" "${POS_ARGS[@]}"
}

set_time_cmd() {
    if /usr/bin/time -v sleep 0 >/dev/null 2>&1; then
        TIME_CMD="/usr/bin/time -v"
    elif gtime -v sleep 0 >/dev/null 2>&1; then
        TIME_CMD="gtime -v"
    else
        echo "error: can't find GNU time" >&2
        exit 1
    fi
}

set_numa_prefix() {
    if [ -n "$NUMA_NODE" ]; then
        if numactl --hardware | grep -q "node $NUMA_NODE "; then
            echo "using numa node: $NUMA_NODE"
            NUMA_PREFIX="numactl --cpunodebind=$NUMA_NODE --membind=$NUMA_NODE "
        else
            echo "error: numa node $NUMA_NODE invalid"
            exit 1
        fi
    else
        NUMA_PREFIX=""
    fi
}

set_tool_vars() {
    local SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

    local TOOL_BIN="$(realpath "$SCRIPT_DIR/../../tools/bin")"
    NAIL="$NUMA_PREFIX $TOOL_BIN/nail"
    MMSEQS="$NUMA_PREFIX $TOOL_BIN/mmseqs"
    BLASTP="$NUMA_PREFIX $TOOL_BIN/blastp"
    PSIBLAST="$NUMA_PREFIX $TOOL_BIN/psiblast"
    MAKEBLASTDB="$NUMA_PREFIX $TOOL_BIN/makeblastdb"
    HMMSEARCH="$NUMA_PREFIX $TOOL_BIN/hmmsearch"
    PHMMER="$NUMA_PREFIX $TOOL_BIN/phmmer"
    DIAMOND="$NUMA_PREFIX $TOOL_BIN/diamond"
    LASTAL="$NUMA_PREFIX $TOOL_BIN/lastal"
    LASTDB="$NUMA_PREFIX $TOOL_BIN/lastdb"
    
    local UTIL="$(realpath "$SCRIPT_DIR/../../util")"
    HMMBALANCE="$NUMA_PREFIX $UTIL/scripts/hmmbalance"
    FASTABALANCE="$NUMA_PREFIX $UTIL/scripts/fastabalance"
}

run_nail() {
    set_time_cmd

    for v in NAIL RESULTS TMP E THREADS QUERY TARGET; do
        check_defined $v $FUNCNAME
    done

    local PREFIX=$1
    local S_ARGS=$2

    local TBL="$RESULTS/$PREFIX.tbl"
    local TIME="$RESULTS/$PREFIX.time"
    local SEEDS="$RESULTS/$PREFIX.seeds"
    local STATS="$RESULTS/$PREFIX.stats"
    local SUMMARY="$RESULTS/$PREFIX.summary"

    echo "running $PREFIX | $S_ARGS"
    echo "   query: $QUERY"
    echo "  target: $TARGET"

    $TIME_CMD -o $TIME \
        $NAIL search \
        -s \
        -t $THREADS \
        --tmp-dir $TMP \
        --stats-results-path $STATS \
        --tbl-out $TBL \
        -E $E \
        $S_ARGS \
        $QUERY $TARGET >> $SUMMARY

    mv $TMP/align_a.tsv $SEEDS
}

run_phmmer() {
    set_time_cmd

    for v in PHMMER RESULTS E THREADS QUERY TARGET; do
        check_defined $v $FUNCNAME
    done

    local PREFIX=$1
    local S_ARGS=$2

    local TBL="$RESULTS/$PREFIX.tbl"
    local DOM="$RESULTS/$PREFIX.domtbl"
    local OUT="$RESULTS/$PREFIX.out"
    local TIME="$RESULTS/$PREFIX.time"

    echo "running $PREFIX | $S_ARGS"
    echo "   query: $QUERY"
    echo "  target: $TARGET"

    $TIME_CMD -o $TIME \
        $PHMMER \
        --cpu $THREADS \
        -E $E \
        $S_ARGS \
        -o /dev/null \
        --tblout $TBL \
        $QUERY $TARGET
}

run_phmmer_split() {
    set_time_cmd

    for v in PHMMER FASTABALANCE RESULTS TMP E THREADS QUERY TARGET; do
        check_defined $v $FUNCNAME
    done

    if (( THREADS % 4 != 0 )); then
        echo "threads: $THREADS"
        echo "threads must be a multiple of 4 (just trust me)"
        exit
    fi

    local PREFIX=$1
    local S_ARGS=$2

    local TBL="$RESULTS/$PREFIX.tbl"
    local DOM="$RESULTS/$PREFIX.domtbl"
    local OUT="$RESULTS/$PREFIX.out"
    local TIME="$RESULTS/$PREFIX.time"

    local N_SPLITS=$(( THREADS / 4 ))
    local SPLIT_THREADS=4
    local SPLIT_DIR=$TMP/query-splits

    $FASTABALANCE $QUERY $N_SPLITS $SPLIT_DIR

    echo "running $PREFIX | $S_ARGS"
    echo "   query: $QUERY"
    echo "  target: $TARGET"

    parallel \
        "${TIME_CMD} -o ${SPLIT_DIR}/{/.}.time \
        $PHMMER \
        --cpu $SPLIT_THREADS \
        -E $E \
        $S_ARGS \
        -o /dev/null \
        --tblout ${SPLIT_DIR}/{/.}.tbl \
        {} ${TARGET}" \
        ::: "${SPLIT_DIR}"/*.fa
    
    cat $SPLIT_DIR/*.tbl > $TBL
    cat $SPLIT_DIR/*.time > $TIME
}

run_hmmsearch() {
    set_time_cmd

    for v in HMMSEARCH RESULTS E THREADS QUERY TARGET; do
        check_defined $v $FUNCNAME
    done

    local PREFIX=$1
    local S_ARGS=$2

    local TBL="$RESULTS/$PREFIX.tbl"
    local DOM="$RESULTS/$PREFIX.domtbl"
    local OUT="$RESULTS/$PREFIX.out"
    local TIME="$RESULTS/$PREFIX.time"

    echo "running $PREFIX | $S_ARGS"
    echo "   query: $QUERY"
    echo "  target: $TARGET"

    $TIME_CMD -o $TIME \
        $HMMSEARCH \
        --cpu $THREADS \
        -E $E \
        $S_ARGS \
        -o /dev/null \
        --tblout $TBL \
        $QUERY $TARGET
}

run_hmmsearch_split() {
    set_time_cmd

    for v in HMMSEARCH HMMBALANCE RESULTS TMP E THREADS QUERY TARGET; do
        check_defined $v $FUNCNAME
    done

    if (( THREADS % 4 != 0 )); then
        echo "threads: $THREADS"
        echo "threads must be a multiple of 4 (just trust me)"
        exit
    fi

    local PREFIX=$1
    local S_ARGS=$2

    local TBL="$RESULTS/$PREFIX.tbl"
    local DOM="$RESULTS/$PREFIX.domtbl"
    local OUT="$RESULTS/$PREFIX.out"
    local TIME="$RESULTS/$PREFIX.time"

    local N_SPLITS=$(( THREADS / 4 ))
    local SPLIT_THREADS=4
    local SPLIT_DIR=$TMP/query-splits

    $HMMBALANCE $QUERY $N_SPLITS $SPLIT_DIR

    echo "running $PREFIX | $S_ARGS"
    echo "   query: $QUERY"
    echo "  target: $TARGET"

    parallel \
        "${TIME_CMD} -o ${SPLIT_DIR}/{/.}.time \
        $HMMSEARCH \
        --cpu $SPLIT_THREADS \
        -E $E \
        $S_ARGS \
        -o /dev/null \
        --tblout ${SPLIT_DIR}/{/.}.tbl \
        {} ${TARGET}" \
        ::: "${SPLIT_DIR}"/*.hmm
    
    cat $SPLIT_DIR/*.tbl > $TBL
    cat $SPLIT_DIR/*.time > $TIME
}

run_mmseqs() {
    set_time_cmd

    for v in MMSEQS RESULTS E THREADS QDB TDB ADB ANNOYING; do
        check_defined $v $FUNCNAME
    done

    local PREFIX=$1
    local S_ARGS=$2

    local TBL="$RESULTS/$PREFIX.tbl"
    local TIME="$RESULTS/$PREFIX.time"

    [ -e $ANNOYING ] && rm -rf $ANNOYING
    [ -e $ADB ] && rm -f $ADB*
    [ -e $ADB.1 ] && rm -f $ADB.*
    
    echo "running $PREFIX | $S_ARGS"
    echo "   query: $QDB"
    echo "  target: $TDB"

    $TIME_CMD -o $TIME \
        $MMSEQS search \
        $QDB $TDB $ADB $ANNOYING \
        --threads $THREADS \
        -e $E \
        $S_ARGS > /dev/null

    $MMSEQS convertalis $QDB $TDB $ADB $TBL --format-mode 0 > /dev/null
}

run_blastp() {
    set_time_cmd

    for v in BLASTP RESULTS E THREADS QUERY_FA TARGET_DB; do
        check_defined $v $FUNCNAME
    done

    local PREFIX=$1
    local S_ARGS=$2

    local TBL="$RESULTS/$PREFIX.tbl"
    local TIME="$RESULTS/$PREFIX.time"

    echo "running $PREFIX | $S_ARGS"
    echo "   query: $QUERY_FA"
    echo "  target: $TARGET_DB"

    $TIME_CMD -o $TIME \
        $BLASTP -query $QUERY_FA \
        -db $TARGET_DB \
        -out $TBL \
        -outfmt 6 \
        -evalue $E \
        -num_threads $THREADS \
        $S_ARGS
}

run_psiblast() {
    set_time_cmd

    for v in PSIBLAST RESULTS E THREADS QUERY_AFA TARGET_DB; do
        check_defined $v $FUNCNAME
    done

    local PREFIX=$1
    local S_ARGS=$2

    local TBL="$RESULTS/$PREFIX.tbl"
    local TIME="$RESULTS/$PREFIX.time"

    echo "running $PREFIX | $S_ARGS"
    echo "   query: $QUERY_AFA"
    echo "  target: $TARGET_DB"

    $TIME_CMD -o $TIME bash -c "\
    for q in $QUERY_AFA/*.afa; do \
      $PSIBLAST -in_msa \$q \
        -db $TARGET_DB \
        -outfmt 6 \
        -evalue $E \
        -num_threads $THREADS \
        -comp_based_stats 1 \
        -num_iterations 1 \
        $S_ARGS >> $TBL; \
    done"
}

run_last() {
    set_time_cmd

    for v in LASTAL RESULTS E THREADS QUERY TARGET_DB; do
        check_defined $v $FUNCNAME
    done

    PREFIX=$1
    S_ARGS=$2

    TBL="$RESULTS/$PREFIX.tbl"
    TIME="$RESULTS/$PREFIX.time"

    echo "running $PREFIX | $S_ARGS"
    echo "   query: $QUERY_FA"
    echo "  target: $TARGET_DB"

    $TIME_CMD -o $TIME \
        $LASTAL $TARGET_DB $QUERY \
        -f BlastTab \
        -P $THREADS \
        $S_ARGS > $TBL

}

run_diamond() {
    set_time_cmd

    for v in DIAMOND RESULTS E THREADS QUERY TARGET_DB; do
        check_defined $v $FUNCNAME
    done

    local PREFIX=$1
    local S_ARGS=$2

    local TBL="$RESULTS/$PREFIX.tbl"
    local TIME="$RESULTS/$PREFIX.time"

    echo "running $PREFIX | $S_ARGS"
    echo "   query: $QUERY"
    echo "  target: $TARGET_DB"

    $TIME_CMD -o $TIME \
        $DIAMOND blastp --query $QUERY \
        --db $TARGET_DB \
        --out $TBL \
        --outfmt 6 \
        --evalue $E \
        --threads $THREADS \
        $S_ARGS > /dev/null 2>&1
}
