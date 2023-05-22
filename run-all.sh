#! /bin/sh

echo hmmer:
time ./scripts/run-hmmsearch.sh
echo

echo mmore: 
time ./scripts/run-mmoreseqs.sh
echo

echo mmeqs:
time ./scripts/run-mmseqs.sh