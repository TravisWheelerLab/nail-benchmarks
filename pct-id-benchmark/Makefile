THREADS     ?= 8
DECOY_RATIO ?= 100

# profmark args
SEED          ?= 67779
N_DECOYS      ?= 5000000
TRAIN_TEST_ID ?= 0.5
MIN_TEST      ?= 10
MAX_TEST      ?= 30

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
TOOLS    := $(realpath $(MAKEFILE_DIR)../tools)
UTIL     := $(realpath $(MAKEFILE_DIR)../util)
DATA     := $(realpath $(MAKEFILE_DIR)../data)
SCRIPTS  := $(MAKEFILE_DIR)/scripts

# hmmer util
PROFMARK   := $(TOOLS)/bin/create-profmark
HMMBUILD   := $(TOOLS)/bin/hmmbuild
HMMBALANCE := $(UTIL)/scripts/hmmbalance

# rust util
RUST     := $(UTIL)/rust
RUST_BIN := $(RUST)/target/release

SPLIT     := $(RUST_BIN)/split-decoys
SPLIT_SRC := $(RUST)/split-decoys/src/main.rs

PID       := $(RUST_BIN)/pid
PID_SRC   := $(RUST)/pid/src

PLOT     := $(RUST_BIN)/pid-plot
PLOT_SRC := $(RUST)/pid-plot/src

# align scripts
NAIL    := $(SCRIPTS)/run-nail.sh
HMMER   := $(SCRIPTS)/run-hmmer.sh
MMSEQS  := $(SCRIPTS)/run-mmseqs.sh
BLAST   := $(SCRIPTS)/run-blast.sh
LAST    := $(SCRIPTS)/run-last.sh
DIAMOND := $(SCRIPTS)/run-diamond.sh

# profmark input 
PM_DIR     := $(MAKEFILE_DIR)/profmark
PM_NAME    := $(PM_DIR)/benchmark
SRC_STO    ?= $(DATA)/pfam.sto
SRC_FA     := $(DATA)/swissprot.fa
SRC_FA_IDX := $(SRC_FA).ssi

# profmark output
TRAIN_MSA := $(PM_NAME).train.msa
TEST_FA   := $(PM_NAME).test.fa
POS_FA    := $(PM_DIR)/positives.fa
DEC_FA    := $(PM_DIR)/decoys.fa
QUERY_STO := $(PM_DIR)/query.sto

PM_OUT    := $(POS_FA) $(DEC_FA) $(QUERY_STO)
PM_REMOVE := $(TEST_FA) $(PM_NAME).neg $(PM_NAME).pos $(PM_NAME).tbl

##########################################
##########################################
##########################################
 
CHECK_PATHS := $(MAKEFILE_DIR) $(TOOLS) $(DATA) $(SRC_FA) $(SRC_STO) $(SPLIT_SRC) $(PID_SRC) $(PLOT_SRC)
.PHONY: $(CHECK_PATHS)
$(CHECK_PATHS):
	@if [ -e "$@" ]; then \
		echo "\"$@\" OK"; \
	else \
		exit 1; \
	fi

check: $(CHECK_PATHS)

$(SRC_FA_IDX):
	@echo indexing \"$(SRC_FA)\"...
	@esl-sfetch --index $(SRC_FA) >/dev/null

$(PM_OUT): $(SRC_FA_IDX)
	@echo running profmark...
	@mkdir -p $(PM_DIR)
	@$(PROFMARK) \
		-S $(SEED) \
		-N $(N_DECOYS) \
		-1 $(TRAIN_TEST_ID) \
		--cluster \
		--mintest $(MIN_TEST) \
		--maxtest $(MAX_TEST) \
		$(PM_NAME) $(SRC_STO) $(SRC_FA)
	
	@$(SPLIT) $(TEST_FA)
	@mv $(TRAIN_MSA) $(QUERY_STO)
	@rm $(PM_REMOVE)
		
profmark: $(PM_OUT)

$(PID): $(PID_SRC)/main.rs
	@cd $(PID_SRC) && cargo build --release

b%: profmark $(PID)
	@mkdir -p b-$*
	@$(PID) b-$* $(POS_FA) $(DEC_FA) $(QUERY_STO) $(SRC_STO) $* $(DECOY_RATIO)
	@$(HMMBUILD) b-$*/query.hmm b-$*/query.sto > /dev/null

run-%:
	@echo "benchmarking $*/ with $(THREADS) threads..."
	@$(NAIL) $* $(THREADS)
	@$(MMSEQS) $* $(THREADS)
	@$(HMMER) $* $(THREADS)
	@$(BLAST) $* $(THREADS)
	@$(LAST) $* $(THREADS)
	@$(DIAMOND) $* $(THREADS)

$(PLOT): $(PLOT_SRC)/main.rs
	@cd $(PLOT_SRC) && cargo build --release

plot-%: $(PLOT)
	$(PLOT) $*

clean:
	@rm -rf tmp/ b-*/results/ benchmark/results/
