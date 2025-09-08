THREADS     ?= 8
DECOY_RATIO ?= 100

# profmark args
SEED          ?= 67779
N_DECOYS      ?= 5000000
TRAIN_TEST_ID ?= 0.5
MIN_TEST      ?= 1
MAX_TEST      ?= 100

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
TRAIN_MSA  := $(PM_NAME).train.msa
TEST_MSA   := $(PM_NAME).test.msa
QUERY_STO  := $(PM_DIR)/query.sto
TARGET_STO := $(PM_DIR)/target.sto
PM_OUT     := $(QUERY_STO) $(TARGET_STO)

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

$(PM_OUT):
	@echo running profmark...
	@mkdir -p $(PM_DIR)
	@$(PROFMARK) \
		-S $(SEED) \
		-N $(N_DECOYS) \
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

pmark:
	@echo running profmark...
	@mkdir -p $(PM_DIR)
	@$(PROFMARK) \
		-S $(SEED) \
		-N $(N_DECOYS) \
		-1 $(TRAIN_TEST_ID) \
		--cluster \
		--onlysplit \
		--mintest $(MIN_TEST) \
		--maxtest $(MAX_TEST) \
		$(PM_NAME) $(SRC_STO)
	
$(PID): $(PID_SRC)/main.rs
	@cd $(PID_SRC) && cargo build --release

b%: profmark $(PID)
	@echo building benchmark...
	@mkdir -p b-$*
	@$(PID) b-$* $(QUERY_STO) $(TARGET_STO) $(SRC_STO) $(SRC_FA) $* $(DECOY_RATIO)
	@echo building query.hmm...
	@$(HMMBUILD) --cpu $(THREADS) b-$*/query.hmm b-$*/query.sto > /dev/null
	@echo done

run-%:
	@echo "benchmarking $*/ with $(THREADS) threads..."
	@$(NAIL)    $* $(THREADS)
	@$(MMSEQS)  $* $(THREADS)
	@$(HMMER)   $* $(THREADS)
	@$(BLAST)   $* $(THREADS)
	@$(LAST)    $* $(THREADS)
	@$(DIAMOND) $* $(THREADS)

$(PLOT): $(PLOT_SRC)/main.rs
	@cd $(PLOT_SRC) && cargo build --release

plot-%: $(PLOT)
	$(PLOT) $*

clean:
	@rm -rf tmp/ b-*/results/ benchmark/results/
