QUERIES=$LONG_SEQ_DIR/query/
TARGETS=$LONG_SEQ_DIR/target/

echo "running nail on long sequence pairs..."
TBL=long-seq.tbl
for ((i=1; i<=6; i++)); do
  Q="$QUERIES${i}.query.fa"
  T="$TARGETS${i}.target.fa"
  nail search --tbl-out tmp.tsv $Q $T
  cat tmp.tsv >> $TBL
  rm tmp.tsv
done
