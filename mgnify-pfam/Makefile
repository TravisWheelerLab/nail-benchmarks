THREADS   ?= 8
NUMA_NODE ?=
SEED      ?= 67779

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(MAKEFILE_DIR)../vars.mk

#############
# BIN SETUP #
#############

PARSE  := $(MGY_PARSE)
SPLIT  := $(FAB)

#############

N_SPLIT := 1000
CUTOFFS := ./cutoffs.txt

# source data
SRC_STO ?= $(DATA)/pfam.sto
SRC_HMM ?= $(DATA)/pfam.hmm 
MGY_FA  ?= $(DATA)/mgnify.fa
MGY_DIR ?= $(DATA)/mgy
MGY_REV_DIR ?= $(DATA)/mgy-rev

#############

small:
	@echo building benchmark...
	$(eval DIR:= benchmark-sm)
	@mkdir -p $(DIR)
	@ln -sf $(SRC_HMM) $(DIR)/query.hmm
	@ln -sf $(SRC_STO) $(DIR)/query.sto
	@$(FASTA_SAMPLE) benchmark-one/target.fa -n 10000 --seed $(SEED) > $(DIR)/target.fa
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
	@ln -sf $(MGY_DIR) $(DIR)/mgy-rev
	@echo done

#############

run-nail-0:
	@./scripts/run-nail.sh \
		--numa_node 0 \
		--threads 48 \
		--a 1 \
		--b 500 \
		--tmp /data/results-mgy/tmp/nail-0/ \
		--results /data/results-mgy/nail/ \
		benchmark-all/

run-nail-1:
	@./scripts/run-nail.sh \
		--numa_node 1 \
		--threads 48 \
		--a 501 \
		--b 1000 \
		--tmp /data/results-mgy/tmp/nail-1/ \
		--results /data/results-mgy/nail/ \
		benchmark-all/

run-rev-0:
	@./scripts/run-rev.sh \
		--numa_node 0 \
		--threads 48 \
		--a 1 \
		--b 500 \
		--tmp /data/results-mgy/tmp/rev-0/ \
		--results /data/results-mgy/rev/ \
		benchmark-all/

run-rev-1:
	@./scripts/run-rev.sh \
		--numa_node 1 \
		--threads 48 \
		--a 501 \
		--b 1000 \
		--tmp /data/results-mgy/tmp/rev-1/ \
		--results /data/results-mgy/rev/ \
		benchmark-all/

run-rev-params:
	@./scripts/params-rev.sh \
		--numa_node 1 \
		--threads 48 \
		--tmp ./tmp-rev/ \
		--results ./results-rev/ \
		benchmark-all/

#############

learn-cutoffs: mgy-parse
	@$(PARSE) cutoffs \
		/data/results-mgy/rev/nail/ \
		/data/results-mgy/rev/mmseqs/ \
		$(SRC_HMM) \
		$(CUTOFFS) \
		-t 48

parse-10: mgy-parse
	@echo "parsing benchmark-all (10)..."
	$(eval DIR:= benchmark-all)
	$(eval RESULTS:= $(DIR)/results/10)
	$(eval ANALYSIS:= $(DIR)/analysis)
	@numactl --cpunodebind=1 --membind=1 $(PARSE) other \
 		--dir \
		-t 10 \
		--print-times \
		$(SRC_HMM) \
		$(CUTOFFS) \
		$(RESULTS)/nail/ \
		$(RESULTS)/mmseqs/ \
		$(RESULTS)/hmmer/ \
		$(ANALYSIS)

parse-100: mgy-parse
	@echo "parsing benchmark-all (100)..."
	$(eval DIR:= benchmark-all)
	$(eval RESULTS:= $(DIR)/results/100)
	$(eval ANALYSIS:= $(DIR)/analysis)
	@numactl --cpunodebind=1 --membind=1 $(PARSE) other \
 		--dir \
		-t 48 \
		--print-times \
		$(SRC_HMM) \
		$(CUTOFFS) \
		$(RESULTS)/nail/ \
		$(RESULTS)/mmseqs/ \
		$(RESULTS)/hmmer/ \
		$(ANALYSIS)

parse-one: mgy-parse
	@echo parsing benchmark-one...
	$(eval DIR:= benchmark-one)
	$(eval RESULTS:= $(DIR)/results)
	$(eval ANALYSIS:= $(DIR)/analysis)
	@mkdir -p $(ANALYSIS)
	@numactl --cpunodebind=1 --membind=1 $(PARSE) other \
		$(SRC_HMM) \
		$(CUTOFFS) \
		$(RESULTS)/nail.prf.tbl \
		$(RESULTS)/mmseqs.prf.tbl \
		$(RESULTS)/hmmer.prf.tbl \
		$(ANALYSIS)
		
#############

run-params-nail: 
	./scripts/params-nail.sh benchmark-one/ --numa_node 1
	
run-params-mmseqs: 
	./scripts/params-mmseqs.sh benchmark-one/ --numa_node 0

params-%: mgy-parse
	$(eval DIR:= $*)
	$(eval P := figures-params/)
	@$(PARSE) params \
		$(DIR)/query.hmm \
		$(DIR)/target.fa \
		$(DIR)/results/hmmer.prf.tbl \
		$(DIR)/results-nail-params \
		$(DIR)/results-mmseqs-params \
		params.txt

