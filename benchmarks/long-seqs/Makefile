QUERIES=data/query/
TARGETS=data/target/
ARGS=--allow-overwrite --f32-p 5
RESULTS=./results
TLEN=$(RESULTS)/target_lengths.txt
QLEN=$(RESULTS)/query_lengths.txt
TBL=$(RESULTS)/long-seq.tbl
POINTS=$(RESULTS)/cells.long.txt

$(RESULTS):
	@mkdir -p $(RESULTS)

$(TLEN): $(RESULTS)
	@rm -rf $(TLEN)
	@esl-seqstat -a $(TARGETS)/1.target.fa | head -n 1 >> $(TLEN)
	@esl-seqstat -a $(TARGETS)/2.target.fa | head -n 1 >> $(TLEN)
	@esl-seqstat -a $(TARGETS)/3.target.fa | head -n 1 >> $(TLEN)
	@esl-seqstat -a $(TARGETS)/4.target.fa | head -n 1 >> $(TLEN)
	@esl-seqstat -a $(TARGETS)/5.target.fa | head -n 1 >> $(TLEN)
	@esl-seqstat -a $(TARGETS)/6.target.fa | head -n 1 >> $(TLEN)

$(QLEN): $(RESULTS)
	@rm -rf $(QLEN) 
	@esl-seqstat -a $(QUERIES)/1.query.fa  | head -n 1 >> $(QLEN)
	@esl-seqstat -a $(QUERIES)/2.query.fa  | head -n 1 >> $(QLEN)
	@esl-seqstat -a $(QUERIES)/3.query.fa  | head -n 1 >> $(QLEN)
	@esl-seqstat -a $(QUERIES)/4.query.fa  | head -n 1 >> $(QLEN)
	@esl-seqstat -a $(QUERIES)/5.query.fa  | head -n 1 >> $(QLEN)
	@esl-seqstat -a $(QUERIES)/6.query.fa  | head -n 1 >> $(QLEN)

$(TBL): $(RESULTS)
	@rm -rf $(TBL)
	@nail search $(QUERIES)/1.query.fa $(TARGETS)/1.target.fa $(ARGS)
	@tail -n 1 results.tbl >> $(TBL)
	@nail search $(QUERIES)/2.query.fa $(TARGETS)/2.target.fa $(ARGS)
	@tail -n 1 results.tbl >> $(TBL)
	@nail search $(QUERIES)/3.query.fa $(TARGETS)/3.target.fa $(ARGS)
	@tail -n 1 results.tbl >> $(TBL)
	@nail search $(QUERIES)/4.query.fa $(TARGETS)/4.target.fa $(ARGS)
	@tail -n 1 results.tbl >> $(TBL)
	@nail search $(QUERIES)/5.query.fa $(TARGETS)/5.target.fa $(ARGS)
	@tail -n 1 results.tbl >> $(TBL)
	@nail search $(QUERIES)/6.query.fa $(TARGETS)/6.target.fa $(ARGS)
	@tail -n 1 results.tbl >> $(TBL)
	@rm results.tbl

parse: $(QLEN) $(TLEN) $(TBL)
	@./parse.sh $(TBL) $(TLEN) $(QLEN) > $(POINTS)
