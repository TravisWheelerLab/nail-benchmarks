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
REV    := $(FAR)

#############

N_SPLIT := 1000
CUTOFFS := cutoffs/

# source data
SRC_STO ?= $(DATA)/pfam.sto
SRC_HMM ?= $(DATA)/pfam.hmm 
MGY_FA  ?= $(DATA)/mgnify.fa
MGY_DIR ?= $(DATA)/mgy
MGY_REV_DIR ?= $(DATA)/mgy-rev

#############
sample: $(FAS)
	@$(FAS) benchmark-one/target.fa -n 10

split: $(FAB)
	@echo splitting $(MGY_FA)...
	@$(SPLIT) -o $(MGY_DIR)/ --n-splits $(N_SPLIT) $(MGY_FA)
	
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
	@ln -sf $(MGY_REV_DIR) $(DIR)/mgy-rev
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

CHECK_REV := numactl --cpunodebind=1 --membind=1 $(PARSE) check-rev \
		--nail-dir   benchmark-all/results/rev/nail/ \
		--mmseqs-dir benchmark-all/results/rev/mmseqs/ \
		--target-dir benchmark-all/mgy/ \
		--query-path $(SRC_HMM) \
		-t $(THREADS) \

check-rev-1: mgy-parse
	@$(CHECK_REV) \
		-n 1 \
		-o $(CUTOFFS)/1

check-rev-10: mgy-parse
	@$(CHECK_REV) \
		-n 10 \
		-o $(CUTOFFS)/10

check-rev-100: mgy-parse
	@$(CHECK_REV) \
		-n 100 \
		-o $(CUTOFFS)/100

check-rev-1000: mgy-parse
	@$(CHECK_REV) \
		-n 1000 \
		-o $(CUTOFFS)/1000

reverse-decoys: far
	@parallel '$(REV) {} -o $(CUTOFFS)/decoys-rev/{/.}.rev.fa' ::: $(CUTOFFS)/decoys/*.fa

decoys:
	@./scripts/run-decoys.sh \
		--numa_node 1 \
		--threads $(THREADS) \
		--tmp $(CUTOFFS)/tmp/ \
		--results $(CUTOFFS)/results/ \
		cutoffs/
		

#############

LEARN_CUTOFFS := numactl --cpunodebind=1 --membind=1 $(PARSE) learn-cutoffs \
		/data/results-mgy/rev/nail/ \
		/data/results-mgy/rev/mmseqs/ \
		$(SRC_HMM) \
		-t $(THREADS) \

learn-cutoffs-1: mgy-parse
	@$(LEARN_CUTOFFS) \
		-n 1 \
		-o $(CUTOFFS)/cutoffs.1.txt

learn-cutoffs-10: mgy-parse
	@$(LEARN_CUTOFFS) \
		-n 10 \
		-o $(CUTOFFS)/cutoffs.10.txt

learn-cutoffs-50: mgy-parse
	@$(LEARN_CUTOFFS) \
		-n 50 \
		-o $(CUTOFFS)/cutoffs.50.txt

learn-cutoffs-100: mgy-parse
	@$(LEARN_CUTOFFS) \
		-n 100 \
		-o $(CUTOFFS)/cutoffs.100.txt

learn-cutoffs-1000: mgy-parse
	@$(LEARN_CUTOFFS) \
		-o $(CUTOFFS)/cutoffs.1000.txt

learn-cutoffs-all:
	@$(MAKE) learn-cutoffs-1    THREADS=$(THREADS)
	@$(MAKE) learn-cutoffs-10   THREADS=$(THREADS)
	@$(MAKE) learn-cutoffs-50   THREADS=$(THREADS)
	@$(MAKE) learn-cutoffs-100  THREADS=$(THREADS)
	@$(MAKE) learn-cutoffs-1000 THREADS=$(THREADS)

#############

RECALL := numactl --cpunodebind=1 --membind=1 \
		  $(PARSE) recall \
		  --query-path $(SRC_HMM) \
		  --cutoffs-path $(CUTOFFS)/cutoffs.1000.txt \
		  --nail-path benchmark-all/results/nail \
		  --mmseqs-path benchmark-all/results/mmseqs \
		  --hmmer-path benchmark-all/results/hmmer \
		  --dir \
		  -t $(THREADS)

recall-1: mgy-parse
	@echo "parsing benchmark-all (0.1%)..."
	@$(RECALL) \
		-n 1 \
		-o benchmark-all/analysis/1
		
recall-10: mgy-parse
	@echo "parsing benchmark-all (1%)..."
	@$(RECALL) \
		-n 10 \
		-o benchmark-all/analysis/10

recall-100: mgy-parse
	@echo "parsing benchmark-all (10%)..."
	@$(RECALL) \
		-n 100 \
		-o benchmark-all/analysis/100

recall-1000: mgy-parse
	@echo "parsing benchmark-all (100%)..."
	@$(RECALL) \
		-n 1000 \
		-o benchmark-all/analysis/1000
		
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

