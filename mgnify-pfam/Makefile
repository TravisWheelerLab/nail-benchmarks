THREADS   ?= 8
NUMA_NODE ?=
SEED      ?= 67779

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(MAKEFILE_DIR)../vars.mk

PARSE_SRC := $(RUST)/mgnify-parse/src/
PARSE     := $(RUST_BIN)/mgnify-parse

SPLIT_SRC := $(RUST)/fab/src/
SPLIT     := $(RUST_BIN)/fab
N_SPLIT   := 1000

# source data
SRC_STO ?= $(DATA)/pfam.sto
SRC_HMM ?= $(DATA)/pfam.hmm 
MGY_FA  ?= $(DATA)/mgnify.fa

# ~1M alignments
small:
	@echo building benchmark...
	$(eval DIR:= benchmark-sm)
	@mkdir -p $(DIR)
	@ln -sf $(SRC_HMM) $(DIR)/query.hmm
	@ln -sf $(SRC_STO) $(DIR)/query.sto
	@$(FASTA_SAMPLE) $(SRC_FA) -n 50 --seed $(SEED) > $(DIR)/target.fa
	@echo done

$(PARSE): $(PARSE_SRC)/main.rs $(BIOIO_SRC)
	@cd $(PARSE_SRC) && cargo build --release

parse-%: $(PARSE)
	@$(PARSE) $*

clean:
	@rm -rf tmp/
