THREADS   ?= 8
NUMA_NODE ?=
SEED      ?= 67779

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(MAKEFILE_DIR)../vars.mk

PARSE_SRC := $(RUST)/bm/src/bin/mgy-parse.rs
PARSE     := $(RUST_BIN)/mgy-parse

PLOT_SRC := $(RUST)/bm/src/bin/mgy-plot.rs
PLOT     := $(RUST_BIN)/mgy-plot

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

run-nail-0:
	@./scripts/run-nail.sh \
		--numa_node 0 \
		--threads 96 \
		--a 1 \
		--b 500 \
		--tmp /data/results-mgy/tmp/nail-0/ \
		--results /data/results-mgy/nail/ \
		benchmark-all/

run-nail-1:
	@./scripts/run-nail.sh \
		--numa_node 1 \
		--threads 96 \
		--a 501 \
		--b 1000 \
		--tmp /data/results-mgy/tmp/nail-1/ \
		--results /data/results-mgy/nail/ \
		benchmark-all/


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

$(PLOT): $(PLOT_SRC) $(BIOIO_SRC)
	@cd $(RUST) && cargo build --release --bin mgy-plot

plot-all: $(PLOT)
	@echo plotting benchmark-all...
	@$(PLOT) -t $(THREADS) -o dist.tbl benchmark-all/results/analysis

run-params-nail: 
	./scripts/params-nail.sh benchmark-one/ --numa_node 1
	
run-params-mmseqs: 
	./scripts/params-mmseqs.sh benchmark-one/ --numa_node 0

params-%: $(PARSE)
	$(eval DIR:= $*)
	$(eval P := figures-params/)
	@$(PARSE) params \
		$(DIR)/query.hmm \
		$(DIR)/target.fa \
		$(DIR)/results/hmmer.prf.tbl \
		$(DIR)/results-nail-params \
		$(DIR)/results-mmseqs-params \
		params.txt

