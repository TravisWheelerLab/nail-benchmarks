THREADS   ?= 8
NUMA_NODE ?=
SEED      ?= 67779

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(MAKEFILE_DIR)../vars.mk

PARSE_SRC := $(RUST)/bm/src/bin/mgy-parse.rs
PARSE     := $(RUST_BIN)/mgy-parse

SPLIT_SRC := $(RUST)/fab/src/
SPLIT     := $(RUST_BIN)/fab
N_SPLIT   := 1000

# source data
SRC_STO ?= $(DATA)/pfam.sto
SRC_HMM ?= $(DATA)/pfam.hmm 
MGY_FA  ?= $(DATA)/mgnify.fa
MGY_DIR ?= $(DATA)/mgy

# ~1M alignments
small:
	@echo building benchmark...
	$(eval DIR:= benchmark-sm)
	@mkdir -p $(DIR)
	@ln -sf $(SRC_HMM) $(DIR)/query.hmm
	@ln -sf $(SRC_STO) $(DIR)/query.sto
	@$(FASTA_SAMPLE) $(SRC_FA) -n 50 --seed $(SEED) > $(DIR)/target.fa
	@echo done

benchmark-one:
	@echo building benchmark...
	$(eval DIR:= benchmark-one)
	@mkdir -p $(DIR)
	@ln -sf $(SRC_HMM) $(DIR)/query.hmm
	@ln -sf $(SRC_STO) $(DIR)/query.sto
	@cp $(MGY_DIR)/1.fa $(DIR)/target.fa
	@echo done

benchmark-all:
	@echo building benchmark...
	$(eval DIR:= benchmark-all)
	@mkdir -p $(DIR)
	@ln -sf $(SRC_HMM) $(DIR)/query.hmm
	@ln -sf $(SRC_STO) $(DIR)/query.sto
	@ln -sf $(MGY_DIR) $(DIR)/mgy
	@echo done

$(PARSE): $(PARSE_SRC) $(BIOIO_SRC)
	@cd $(RUST) && cargo build --release --bin mgy-parse

parse-all: $(PARSE)
	@echo parsing benchmark-all...
	$(eval DIR:= benchmark-all)
	$(eval RESULTS:= $(DIR)/results)
	$(eval HMMER:= $(RESULTS)/hmmer)
	$(eval NAIL:= $(RESULTS)/nail)
	$(eval MMSEQS:= $(RESULTS)/mmseqs)
	$(eval ANALYSIS:= $(RESULTS)/analysis)
	@mkdir -p $(ANALYSIS)
	@seq 1 1000 | parallel -j$(THREADS) $(PARSE) \
		$(DIR)/query.hmm \
		$(DIR)/mgy/{}.fa \
		$(HMMER)/hmmer.{}.prf.tbl \
		$(HMMER)/hmmer.{}.prf.domtbl \
		$(NAIL)/nail.{}.prf.tbl \
		$(NAIL)/nail.{}.prf.stats \
		$(MMSEQS)/mmseqs.{}.prf.tbl \
		$(ANALYSIS)/ga_hits.{}.tbl


parse-1: $(PARSE)
	@echo parsing benchmark-all...
	$(eval DIR:= benchmark-all)
	$(eval RESULTS:= $(DIR)/results)
	$(eval HMMER:= $(RESULTS)/hmmer)
	$(eval NAIL:= $(RESULTS)/nail)
	$(eval MMSEQS:= $(RESULTS)/mmseqs)
	$(eval ANALYSIS:= $(RESULTS)/analysis)
	@mkdir -p $(ANALYSIS)
	$(PARSE) \
		$(DIR)/query.hmm \
		$(DIR)/mgy/1.fa \
		$(HMMER)/hmmer.1.prf.tbl \
		$(HMMER)/hmmer.1.prf.domtbl \
		$(NAIL)/nail.1.prf.tbl \
		$(NAIL)/nail.1.prf.stats \
		$(MMSEQS)/mmseqs.1.prf.tbl \
		$(ANALYSIS)/ga_hits.1.tbl

clean:
	@rm -rf tmp/
