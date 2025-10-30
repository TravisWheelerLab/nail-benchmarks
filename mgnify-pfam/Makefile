THREADS   ?= 8
NUMA_NODE ?=
SEED      ?= 67779

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
TOOLS    := $(realpath $(MAKEFILE_DIR)../tools)
UTIL     := $(realpath $(MAKEFILE_DIR)../util)
DATA     := $(realpath $(MAKEFILE_DIR)../data)
SCRIPTS  := $(MAKEFILE_DIR)/scripts

FASTASAMPLE := $(UTIL)/scripts/fastasample.py

# source data
SRC_STO ?= $(DATA)/pfam.sto
# pfam = 20,795 profiles
SRC_HMM ?= $(DATA)/pfam.hmm 
# mgy_proteins 1 = 100M sequences
# SRC_FA  ?= $(DATA)/mgy_proteins_1.fa
SRC_FA  ?= $(DATA)/1m.mgy.fa

# ~1M alignments
small:
	@echo building benchmark...
	$(eval DIR:= benchmark-sm)
	@mkdir -p $(DIR)
	@ln -sf $(SRC_HMM) $(DIR)/query.hmm
	@ln -sf $(SRC_STO) $(DIR)/query.sto
	@$(FASTASAMPLE) $(SRC_FA) -n 50 --seed $(SEED) > $(DIR)/target.fa
	@echo done

# ~10M alignments
medium:
	@echo building benchmark...
	$(eval DIR:= benchmark-med)
	@mkdir -p $(DIR)
	@ln -sf $(SRC_HMM) $(DIR)/query.hmm
	@ln -sf $(SRC_STO) $(DIR)/query.sto
	@$(FASTASAMPLE) $(SRC_FA) -n 500 --seed $(SEED) > $(DIR)/target.fa
	@echo done

# ~100M alignments
large:
	@echo building benchmark...
	$(eval DIR:= benchmark-lg)
	@mkdir -p $(DIR)
	@ln -sf $(SRC_HMM) $(DIR)/query.hmm
	@ln -sf $(SRC_STO) $(DIR)/query.sto
	@$(FASTASAMPLE) $(SRC_FA) -n 5000 --seed $(SEED) > $(DIR)/target.fa
	@echo done

benchmark:
	@echo building benchmark...
	@echo done

run-%:
	@echo "benchmarking $*/ with $(THREADS) threads..."

clean:
	@rm -rf tmp/ benchmark*/results*/
