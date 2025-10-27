THREADS   ?= 8
NUMA_NODE ?=
# profmark args
SEED          ?= 67779
TRAIN_TEST_ID ?= 0.30
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
HMMEMIT    := $(TOOLS)/bin/hmmemit
HMMBALANCE := $(UTIL)/scripts/hmmbalance

# rust util
RUST     := $(UTIL)/rust
RUST_BIN := $(RUST)/target/release

# benchmark creation
PID       := $(RUST_BIN)/pid
PID_SRC   := $(RUST)/pid/src

# plotting
FIGURES   := $(MAKEFILE_DIR)/figures

PARSE     := $(RUST_BIN)/pid-plot
PARSE_SRC := $(RUST)/pid-plot/src

PLOT_PID  := $(SCRIPTS)/plot-pid.py
PLOT_ROC  := $(SCRIPTS)/plot-roc.py
PLOT_TIME := $(SCRIPTS)/plot-time.py

DATA_PID  := $(FIGURES)/pid.txt
DATA_ROC  := $(FIGURES)/roc.txt
DATA_TIME := $(FIGURES)/time.txt

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
 
CHECK_PATHS := $(MAKEFILE_DIR) $(TOOLS) $(DATA) $(SRC_FA) $(SRC_STO) $(SPLIT_SRC) $(PID_SRC) $(PARSE_SRC)
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
	$(PROFMARK) \
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

small: profmark $(PID)
	@echo building benchmark...
	@mkdir -p benchmark-sm/
	@$(PID) benchmark-sm/ $(QUERY_STO) $(TARGET_STO) $(SRC_STO) $(SRC_FA) 100
	@echo building query.hmm...
	@$(HMMBUILD) --cpu $(THREADS) benchmark-sm/query.hmm benchmark-sm/query.sto > /dev/null
	@$(HMMEMIT) -c benchmark-sm/query.hmm > benchmark-sm/query.cons.fa
	@echo done

medium: profmark $(PID)
	@echo building benchmark...
	@mkdir -p benchmark-med/
	@$(PID) benchmark-med/ $(QUERY_STO) $(TARGET_STO) $(SRC_STO) $(SRC_FA) 500
	@echo building query.hmm...
	@$(HMMBUILD) --cpu $(THREADS) benchmark-med/query.hmm benchmark-med/query.sto > /dev/null
	@$(HMMEMIT) -c benchmark-med/query.hmm > benchmark-med/query.cons.fa
	@echo done

large: profmark $(PID)
	@echo building benchmark...
	@mkdir -p benchmark-lg/
	@$(PID) benchmark-lg/ $(QUERY_STO) $(TARGET_STO) $(SRC_STO) $(SRC_FA) 1000
	@echo building query.hmm...
	@$(HMMBUILD) --cpu $(THREADS) benchmark-lg/query.hmm benchmark-lg/query.sto > /dev/null
	@$(HMMEMIT) -c benchmark-lg/query.hmm > benchmark-lg/query.cons.fa
	@echo done

benchmark: profmark $(PID)
	@echo building benchmark...
	@mkdir -p benchmark/
	@$(PID) benchmark/ $(QUERY_STO) $(TARGET_STO) $(SRC_STO) $(SRC_FA)
	@echo building query.hmm...
	@$(HMMBUILD) --cpu $(THREADS) benchmark/query.hmm benchmark/query.sto > /dev/null
	@$(HMMEMIT) -c benchmark/query.hmm > benchmark/query.cons.fa
	@echo done

run-%:
	@echo "benchmarking $*/ with $(THREADS) threads..."
	@$(HMMER)   $* $(THREADS) $(NUMA_NODE) 
	@$(NAIL)    $* $(THREADS) $(NUMA_NODE)
	@$(MMSEQS)  $* $(THREADS) $(NUMA_NODE)
	@$(LAST)    $* $(THREADS) $(NUMA_NODE)
	@$(BLAST)   $* $(THREADS) $(NUMA_NODE)
	@$(DIAMOND) $* $(THREADS) $(NUMA_NODE)

$(PARSE): $(PARSE_SRC)/main.rs
	@cd $(PARSE_SRC) && cargo build --release

plot-%: $(PARSE)
	@$(PARSE) $*
	@$(PLOT_PID) $(DATA_PID)
	@$(PLOT_ROC) $(DATA_ROC)
	@$(PLOT_TIME) $(DATA_TIME)

just-plot-%: $(PARSE)
	@$(PLOT_PID) $(DATA_PID)
	@$(PLOT_ROC) $(DATA_ROC)
	@$(PLOT_TIME) $(DATA_TIME)

clean:
	@rm -rf tmp/ benchmark*/results*/ benchmark*/p7-queryDB
