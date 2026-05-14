#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/init.sh"

set_default THREADS 8
set_default E 10

check_defined TOOL_BIN
check_defined BM_DIR 
check_defined TMP
check_defined RESULTS

mkdir -p $RESULTS/nail
mkdir -p $RESULTS/mmseqs
mkdir -p $RESULTS/hmmer

for ((i=1; i<=${THREADS}; i++)); do
    mkdir -p $TMP/nail/${i}
    mkdir -p $TMP/mmseqs/${i}
done

HMM=$BM_DIR/queries-hmm
STO=$BM_DIR/queries-sto
DECOYS=$BM_DIR/decoys
DECOYS_REV=$BM_DIR/decoys-rev

export NAIL="$TOOL_BIN/nail"
export MMSEQS="$TOOL_BIN/mmseqs"
export HMMSEARCH="$TOOL_BIN/hmmsearch"

export TMP
export HMM
export STO
export DECOYS
export DECOYS_REV
export RESULTS

# # ---

# export ARGS="-t 1 --allow-overwrite --mmseqs-s 12.0 --mmseqs-max-seqs 1000000000"
# time ${NUMA_PREFIX} parallel -j "${THREADS}" '
#     Q="$HMM/{/.}.hmm"
#     T="$DECOYS/{/.}.fa"
#     T_R="$DECOYS_REV/{/.}.rev.fa"
#     TBL="${RESULTS}/nail/{/.}.tbl"
#     TBL_R="${RESULTS}/nail/{/.}.rev.tbl"
#     D="${TMP}/nail/{%}" 

#     $NAIL search \
#         $Q $T \
#         $ARGS \
#         --tmp-dir $D \
#         --tbl-out $TBL \
#         > /dev/null

#     $NAIL search \
#         $Q $T_R \
#         $ARGS \
#         --tmp-dir $D \
#         --tbl-out $TBL_R \
#         > /dev/null
# ' ::: "${DECOYS}"/*

# # ---

# export ARGS="--threads 1 -s 12.0 --max-seqs 1000000000"
# time ${NUMA_PREFIX} parallel -j "${THREADS}" '
#     Q="$STO/{/.}.sto"
#     T="$DECOYS/{/.}.fa"
#     T_R="$DECOYS_REV/{/.}.rev.fa"
#     TBL="${RESULTS}/mmseqs/{/.}.tbl"
#     TBL_R="${RESULTS}/mmseqs/{/.}.rev.tbl"
#     D="${TMP}/mmseqs/{%}" 

#     [ -e $D ] && rm -rf $D
#     mkdir $D

#     ANNOYING="$D/annoying"
    
#     MDB=$D/msaDB
#     $MMSEQS convertmsa "$Q" $MDB --identifier-field 0 > /dev/null
    
#     QDB=$D/queryDB
#     $MMSEQS msa2profile $MDB $QDB --match-mode 1 > /dev/null

#     TDB=$D/targetDB
#     $MMSEQS createdb "$T" $TDB > /dev/null

#     ADB=$D/alignDB

#     $MMSEQS search \
#         $QDB $TDB $ADB $ANNOYING \
#         $ARGS \
#         > /dev/null

#     $MMSEQS convertalis $QDB $TDB $ADB $TBL --format-mode 0 > /dev/null

#     rm $TDB*
#     rm $ADB*
#     rm -rf $ANNOYING
#     $MMSEQS createdb "$T_R" $TDB > /dev/null

#     $MMSEQS search \
#         $QDB $TDB $ADB $ANNOYING \
#         $ARGS \
#         > /dev/null

#     $MMSEQS convertalis $QDB $TDB $ADB $TBL_R --format-mode 0 > /dev/null
#  ' ::: "${DECOYS}"/*


# ---

export ARGS="--cpu 1"
time ${NUMA_PREFIX} parallel -j "${THREADS}" '
    Q="$HMM/{/.}.hmm"
    T="$DECOYS/{/.}.fa"
    T_R="$DECOYS_REV/{/.}.rev.fa"
    TBL="${RESULTS}/hmmer/{/.}.tbl"
    TBL_R="${RESULTS}/hmmer/{/.}.rev.tbl"
    D="${TMP}/hmmer/{%}" 

    $HMMSEARCH \
        $ARGS \
        --tblout $TBL \
        $Q $T \
        > /dev/null

    $HMMSEARCH \
        $ARGS \
        --tblout $TBL_R \
        $Q $T_R \
        > /dev/null
 ' ::: "${DECOYS}"/*
