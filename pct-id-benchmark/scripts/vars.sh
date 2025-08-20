set_vars() {
    if [ "$#" == 0 ]; then
        echo "usage: ./run-<tool>.sh <benchmark-dir> [threads]"
        exit
    fi
    
    if [ -n "$2" ]; then
        export THREADS=$2
    else
        export THREADS=8
    fi

    export DIR=$1
    export QUERY_HMM=$DIR/query.hmm
    export QUERY_MSA=$DIR/query.sto
    export QUERY_FA=$DIR/query.fa
    export QUERY_CONS_FA=$DIR/query.cons.fa
    export TARGET=$DIR/target.fa

    export VAR E=1e9

    export RESULTS=$DIR/results/
    mkdir -p $RESULTS
}
