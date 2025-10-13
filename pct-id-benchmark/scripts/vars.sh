set_vars() {
    if [ "$#" == 0 ]; then
        echo "usage: ./run-<tool>.sh <benchmark-dir> [threads] [numa node]"
        exit
    fi
    
    if [ -n "$2" ]; then
        export THREADS=$2
    else
        export THREADS=8
    fi

    if [ -n "$3" ]; then
        if numactl --hardware | grep -q "node $3 "; then
            echo "using numa node: $3"
            export NUMA_PREFIX="numactl --cpunodebind=$3 --membind=$3 "
        else
            echo "error: numa node $3 invalid"
            exit 1
        fi
    else
        export NUMA_PREFIX=""
    fi

    export DIR=$1
    export BENCH_TBL=$DIR/benchmark.tbl
    export QUERY_HMM=$DIR/query.hmm
    export QUERY_MSA=$DIR/query.sto
    export QUERY_FA=$DIR/query.fa
    export QUERY_CONS_FA=$DIR/query.cons.fa
    export QUERY_AFA=$DIR/afa
    export TARGET=$DIR/target.fa
    export VAR E=1000
    export RESULTS=$DIR/results/
    mkdir -p $RESULTS
}
