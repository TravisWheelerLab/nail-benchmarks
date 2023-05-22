#! /bin/sh

NAME=benchmark
POS=./$NAME/$NAME.pos

SCRIPT=./scripts/recall.py

OUTPUT_DIR=./output

HMMER_RESULTS=$OUTPUT_DIR/hmmer/results.sorted
HMMER_TARGET_COL=0
HMMER_QUERY_COL=2

MMORE_RESULTS=$OUTPUT_DIR/mmoreseqs/results.sorted
MMORE_TARGET_COL=0
MMORE_QUERY_COL=1

NUM_TRUE_POSITIVES="$(wc -l < $POS)"

echo hmmer: $(python3 $SCRIPT $HMMER_RESULTS $HMMER_TARGET_COL $HMMER_QUERY_COL $NUM_TRUE_POSITIVES)
echo mmore: $(python3 $SCRIPT $MMORE_RESULTS $MMORE_TARGET_COL $MMORE_QUERY_COL $NUM_TRUE_POSITIVES)