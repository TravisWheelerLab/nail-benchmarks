MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
THREADS ?= 8

TOOLS    := $(realpath $(MAKEFILE_DIR)../tools)
UTIL     := $(realpath $(MAKEFILE_DIR)../util)
RUST     := $(UTIL)/rust
RUST_BIN := $(RUST)/target/release
DATA     := $(realpath $(MAKEFILE_DIR)../data)
SCRIPTS  := $(MAKEFILE_DIR)/scripts

PROFMARK := $(TOOLS)/bin/create-profmark
HMMBUILD := $(TOOLS)/bin/hmmbuild
SPLIT    := $(RUST_BIN)/split-decoys
PLOT     := $(RUST_BIN)/plot-pid

NAIL   := $(SCRIPTS)/run-nail.sh
HMMER  := $(SCRIPTS)/run-hmmer.sh
MMSEQS := $(SCRIPTS)/run-mmseqs.sh
BLAST  := $(SCRIPTS)/run-blast.sh
LAST   := $(SCRIPTS)/run-last.sh

PM_DIR     := $(MAKEFILE_DIR)/profmark
PM_NAME    := $(PM_DIR)/benchmark
SRC_STO    ?= $(DATA)/pfam.sto
SRC_FA     := $(DATA)/swissprot.fa
SRC_FA_IDX := $(SRC_FA).ssi

# profmark args
SEED          ?= 67779
N_DECOYS      ?= 5000000
TRAIN_TEST_ID ?= 0.5
MIN_TEST      ?= 10
MAX_TEST      ?= 30

# profmark files
TRAIN_MSA := $(PM_NAME).train.msa
TEST_FA   := $(PM_NAME).test.fa
POS_FA    := $(PM_DIR)/positives.fa
DEC_FA    := $(PM_DIR)/decoys.fa
QUERY_STO := $(PM_DIR)/query.sto

PM_OUT    := $(POS_FA) $(DEC_FA) $(QUERY_STO)
PM_REMOVE := $(TEST_FA) $(PM_NAME).neg $(PM_NAME).pos $(PM_NAME).tbl

PID         := $(RUST_BIN)/pid
PID_SRC     := $(RUST)/pid/src
PID_ARGS    := $(POS_FA) $(DEC_FA) $(QUERY_STO) $(SRC_STO)
DECOY_RATIO ?= 100

CHECK_PATHS := $(MAKEFILE_DIR) $(TOOLS) $(DATA) $(SRC_FA)	

.PHONY: $(CHECK_PATHS)
$(CHECK_PATHS):
	@if [ ! -e "$@" ]; then \
		echo "Error: \"$@\" does not exist"; \
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
	@$(PID) b-$* $(PID_ARGS) $* $(DECOY_RATIO)
	@$(HMMBUILD) b-$*/query.hmm b-$*/query.sto > /dev/null

run-%:
	@echo "benchmarking $*/ with $(THREADS) threads..."
	@$(NAIL) $* $(THREADS)
	@$(MMSEQS) $* $(THREADS)
	@$(HMMER) $* $(THREADS)

plot-%:
	$(PLOT) $*

clean:
	@rm -rf tmp-nail/ tmp-mmseqs/ b-*/results/ benchmark/results/
