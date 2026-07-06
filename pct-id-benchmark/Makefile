MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(MAKEFILE_DIR)../vars.mk

# configurable args
THREADS   ?= 8
NUMA_NODE ?=
FIGURES   ?= figures/

# profmark args
SEED          ?= 67779
TRAIN_TEST_ID ?= 0.30
MIN_TEST      ?= 10
MAX_TEST      ?= 30

#############
# BIN SETUP #
#############

CREATE := $(PID)
PARSE  := $(PID_PARSE)

#############

SCRIPTS  := $(MAKEFILE_DIR)/scripts

# plotting
PLOT_PID    := $(SCRIPTS)/plot_pid.py
PLOT_ROC    := $(SCRIPTS)/plot_roc.py
PLOT_TIME   := $(SCRIPTS)/plot_time.py
PLOT_CELLS  := $(SCRIPTS)/plot_cell_frac.py
PLOT_SCORE  := $(SCRIPTS)/plot_score.py
PLOT_PARAMS := $(SCRIPTS)/plot_params.py

FIGURE_DATA :=$(FIGURES)/00data/

OUT_PID    := $(FIGURES)/pid.pdf
OUT_ROC    := $(FIGURES)/roc.pdf
OUT_TIME   := $(FIGURES)/time.pdf
OUT_CELLS  := $(FIGURES)/cells.pdf
OUT_SCORE  := $(FIGURES)/score.pdf
OUT_PARAMS := $(FIGURES)/params.pdf

# align scripts
NAIL    := $(SCRIPTS)/run-nail.sh
HMMER   := $(SCRIPTS)/run-hmmer.sh
MMSEQS  := $(SCRIPTS)/run-mmseqs.sh
LAST    := $(SCRIPTS)/run-last.sh
BLAST   := $(SCRIPTS)/run-blast.sh
DIAMOND := $(SCRIPTS)/run-diamond.sh

# profmark input 
PM_DIR     := $(MAKEFILE_DIR)/profmark
PM_NAME    := $(PM_DIR)/benchmark
SRC_STO    ?= $(DATA)/pfam.sto
SRC_FA     := $(DATA)/swissprot.fa
SRC_FA_IDX := $(SRC_FA).ssi

# profmark output
TRAIN_MSA  := $(PM_NAME).train.msa
TEST_MSA   := $(PM_NAME).test.msa
QUERY_STO  := $(PM_DIR)/query.sto
TARGET_STO := $(PM_DIR)/target.sto
PM_OUT     := $(QUERY_STO) $(TARGET_STO)

##########################################
##########################################
##########################################

$(PM_OUT):
	@echo running profmark...
	@mkdir -p $(PM_DIR)
	@$(PROFMARK) \
		-S $(SEED) \
		-1 $(TRAIN_TEST_ID) \
		--cluster \
		--onlysplit \
		--mintest $(MIN_TEST) \
		--maxtest $(MAX_TEST) \
		$(PM_NAME) $(SRC_STO)
	
	@mv $(TRAIN_MSA) $(QUERY_STO)
	@mv $(TEST_MSA) $(TARGET_STO)
	@rm $(PM_NAME).tbl
		
profmark: $(PM_OUT)

small: profmark $(CREATE)
	@echo building benchmark...
	@mkdir -p benchmark-sm/
	@$(CREATE) benchmark-sm/ $(QUERY_STO) $(TARGET_STO) $(SRC_STO) $(SRC_FA) 100
	@echo building query.hmm...
	@$(HMMBUILD) --cpu $(THREADS) benchmark-sm/query.hmm benchmark-sm/query.sto > /dev/null
	@$(HMMEMIT) -c benchmark-sm/query.hmm > benchmark-sm/query.cons.fa
	@echo done

medium: profmark $(CREATE)
	@echo building benchmark...
	@mkdir -p benchmark-med/
	@$(CREATE) benchmark-med/ $(QUERY_STO) $(TARGET_STO) $(SRC_STO) $(SRC_FA) 500
	@echo building query.hmm...
	@$(HMMBUILD) --cpu $(THREADS) benchmark-med/query.hmm benchmark-med/query.sto > /dev/null
	@$(HMMEMIT) -c benchmark-med/query.hmm > benchmark-med/query.cons.fa
	@echo done

large: profmark $(CREATE)
	@echo building benchmark...
	@mkdir -p benchmark-lg/
	@$(CREATE) benchmark-lg/ $(QUERY_STO) $(TARGET_STO) $(SRC_STO) $(SRC_FA) 1000
	@echo building query.hmm...
	@$(HMMBUILD) --cpu $(THREADS) benchmark-lg/query.hmm benchmark-lg/query.sto > /dev/null
	@$(HMMEMIT) -c benchmark-lg/query.hmm > benchmark-lg/query.cons.fa
	@echo done

benchmark: profmark $(CREATE)
	@echo building benchmark...
	@mkdir -p benchmark/
	@$(CREATE) benchmark/ $(QUERY_STO) $(TARGET_STO) $(SRC_STO) $(SRC_FA)
	@echo building query.hmm...
	@$(HMMBUILD) --cpu $(THREADS) benchmark/query.hmm benchmark/query.sto > /dev/null
	@$(HMMEMIT) -c benchmark/query.hmm > benchmark/query.cons.fa
	@echo done

run-%:
	@echo "benchmarking $*/ with $(THREADS) threads..."
	$(eval RESULTS := $*/results/)
	$(eval N := $(if $(NUMA_NODE),--numa_node $(NUMA_NODE),))
	$(eval ARGS := --threads $(THREADS) --results $(RESULTS) $(N) $*)
	@rm -rf   $(RESULTS)
	@mkdir -p $(RESULTS)
	@$(HMMER)   $(ARGS) 
	@$(NAIL)    $(ARGS) 
	@$(MMSEQS)  $(ARGS) 
	@$(LAST)    $(ARGS) 
	@$(DIAMOND) $(ARGS) 
	@$(BLAST)   $(ARGS)

run-sparse-%:
	@echo "benchmarking $*/ with $(THREADS) threads..."
	$(eval N := $(if $(NUMA_NODE),--numa_node $(NUMA_NODE),))
	$(eval ARGS := --threads $(THREADS) --results $*/results-sparse/ $(N) $*)
	@mkdir -p $*/results-sparse/
	@$(SCRIPTS)/run-sparse.sh $(ARGS) 

plot-recall-%: $(PARSE)
	@$(PARSE) recall $*/benchmark.tbl $*/results/ $(FIGURE_DATA)
	@$(PLOT_PID)  $(FIGURE_DATA)/pid.txt  $(OUT_PID)
	@$(PLOT_ROC)  $(FIGURE_DATA)/roc.txt  $(OUT_ROC)
	@$(PLOT_TIME) $(FIGURE_DATA)/time.txt $(OUT_TIME)

plot-score-%: $(PARSE)
	@$(PARSE) score $*/results-full-sparse/nail-s9.0-full.prf.tbl $*/results-full-sparse/nail-s9.0-sparse.prf.tbl $(FIGURE_DATA)
	# @$(PARSE) score $*/results-sparse/nail-s12.0-full.prf.tbl $*/results/nail-s12.0-ms2000.prf.tbl $(FIGURE_DATA)
	# @$(PARSE) score $*/results-sparse/nail-s12.0-full.prf.tbl $*/results-sparse/nail-s12.0-cells.prf.tbl $(FIGURE_DATA)
	@$(PLOT_SCORE) $(FIGURE_DATA)/score.txt $(OUT_SCORE)

plot-cells-%: $(PARSE)
	$(eval LONG_SEQ:= $(if $(wildcard ../long-seqs/results/cells.long.txt),--long_hits ../long-seqs/results/cells.long.txt,))
	@$(PARSE) cells $*/results-sparse/nail-s12.0-cells.prf.tbl $*/query.hmm $*/target.fa $(FIGURE_DATA)
	@$(PLOT_CELLS) $(FIGURE_DATA)/cells.true.txt $(FIGURE_DATA)/cells.decoy.txt $(OUT_CELLS) $(LONG_SEQ)

plot-all-%:
	@$(MAKE) plot-recall-$*
	@$(MAKE) plot-score-$*
	@$(MAKE) plot-cells-$*

table-score: $(PARSE)
	@$(PARSE) table\
		--benchmark-tbl ./benchmark/benchmark.tbl \
		--results-dir ./benchmark/results/ \
		--out-path results.score.tbl

table-e: $(PARSE)
	@$(PARSE) table\
		--benchmark-tbl ./benchmark/benchmark.tbl \
		--results-dir ./benchmark/results/ \
		--out-path results.evalue.tbl \
		--e-value


params-%: $(PARSE)
	$(eval P := figures-params/)
	$(eval OVERLAY := $(if $(wildcard figures/time.txt),--overlay figures/time.txt,))
	@mkdir -p $(P)/tmp/
	@$(PARSE) $*/benchmark.tbl $*/results-nail-params/ $(P)/tmp/
	@mv $(P)/tmp/time.txt $(P)/nail.txt
	@$(PARSE) $*/benchmark.tbl $*/results-mmseqs-params/ $(P)/tmp/
	@mv $(P)/tmp/time.txt $(P)/mmseqs.txt
	@rm -rf $(P)/tmp/
	@$(PLOT_PARAMS) $(P)/nail.txt $(P)/mmseqs.txt $(P)/params.pdf
	@$(PLOT_PARAMS) $(OVERLAY) $(P)/nail.txt $(P)/mmseqs.txt $(P)/overlay.pdf

clean:
	@rm -rf tmp/ benchmark*/results*/ benchmark*/p7-queryDB
