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

run-mmseqs-single:
	@./scripts/run-mmseqs.sh \
		--numa_node 1 \
		--threads 48 \
		--a 1 \
		--b 1 \
		--tmp /data/results-mgy/tmp/mmseqs-single/ \
		--results /data/results-mgy/mmseqs-single/ \
		benchmark-all/

run-mmseqs-s7.5-a:
	@./scripts/run-mmseqs.sh \
		--numa_node 0 \
		--threads 48 \
		--s 7.5 \
		--a 1 \
		--b 500 \
		--tmp /data/results-mgy/tmp/mmseqs-a/ \
		--results /data/results-mgy/mmseqs-s7.5/ \
		benchmark-all/

run-mmseqs-s7.5-b:
	@./scripts/run-mmseqs.sh \
		--numa_node 1 \
		--threads 48 \
		--s 7.5 \
		--a 501 \
		--b 1000 \
		--tmp /data/results-mgy/tmp/mmseqs-b/ \
		--results /data/results-mgy/mmseqs-s7.5/ \
		benchmark-all/

run-hmmer-single:
	@./scripts/run-hmmer.sh \
		--numa_node 1 \
		--threads 48 \
		--a 1 \
		--b 1 \
		--tmp /data/results-mgy/tmp/hmmer-single/ \
		--results /data/results-mgy/hmmer-single/ \
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

runtime-extrap: 
	@./scripts/runtime-extrap.sh \
		--numa_node $(NUMA_NODE) \
		--threads $(THREADS) \
		--tmp benchmark-all/tmp-runtime-extrap/ \
		benchmark-all/

full-sparse: 
	@./scripts/full-sparse.sh \
		--numa_node $(NUMA_NODE) \
		--threads $(THREADS) \
		--tmp benchmark-all/tmp-full-sparse/ \
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
		--numa_node $(NUMA_NODE) \
		--threads $(THREADS) \
		--tmp $(CUTOFFS)/tmp/ \
		--results $(CUTOFFS)/results/ \
		cutoffs/
		

#############

LEARN_CUTOFFS := numactl --cpunodebind=1 --membind=1 $(PARSE) learn-cutoffs \
		$(CUTOFFS)/results/nail/ \
		$(CUTOFFS)/results/mmseqs/ \
		$(CUTOFFS)/results/hmmer/ \
		$(SRC_HMM) \
		-t $(THREADS) \

learn-cutoffs: mgy-parse
	@$(LEARN_CUTOFFS) \
		-o $(CUTOFFS)/cutoffs.txt

#############

RECALL := numactl --cpunodebind=1 --membind=1 \
		  $(PARSE) recall \
		  --query-path $(SRC_HMM) \
		  --cutoffs-path $(CUTOFFS)/cutoffs.1000.txt \
		  --hmmer-path benchmark-all/results/hmmer \
		  --dir \
		  -t $(THREADS)
		
recall-n9.0-m7.5: mgy-parse
	@echo "parsing benchmark-all (100%)..."
	@$(RECALL) \
		-n 1000 \
		--nail-path benchmark-all/results/nail-s9.0 \
		--mmseqs-path benchmark-all/results/mmseqs-s7.5 \
		-o benchmark-all/analysis/n9.0-m7.5

recall-n10.0-m7.5: mgy-parse
	@echo "parsing benchmark-all (100%)..."
	@$(RECALL) \
		-n 1000 \
		--nail-path benchmark-all/results/nail-s10.0 \
		--mmseqs-path benchmark-all/results/mmseqs-s7.5 \
		-o benchmark-all/analysis/n10.0-m7.5

recall-n12.0-m7.5: mgy-parse
	@echo "parsing benchmark-all (100%)..."
	@$(RECALL) \
		-n 1000 \
		--nail-path benchmark-all/results/nail-s12.0 \
		--mmseqs-path benchmark-all/results/mmseqs-s7.5 \
		-o benchmark-all/analysis/n12.0-m7.5

recall-n12.0-m12.0: mgy-parse
	@echo "parsing benchmark-all (100%)..."
	@$(RECALL) \
		-n 1000 \
		--nail-path benchmark-all/results/nail-s12.0 \
		--mmseqs-path benchmark-all/results/mmseqs-s12.0 \
		-o benchmark-all/analysis/n12.0-m12.0
		
#############

run-params-nail: 
	./scripts/params-nail.sh benchmark-one/ --numa_node 1
	
run-params-mmseqs: 
	./scripts/params-mmseqs.sh benchmark-one/ --numa_node 0

params: mgy-parse
	@$(PARSE) params \
		--query-path   benchmark-one/query.hmm \
		--cutoffs-path $(CUTOFFS)/cutoffs.1000.txt \
		--hmmer-path   benchmark-one/results/hmmer.prf.tbl \
		--hmmer-dom-path   benchmark-one/results/hmmer.prf.domtbl \
		--nail-path    benchmark-one/results-nail-params \
		--mmseqs-path  benchmark-one/results-mmseqs-params \
		--out-path     params.txt



#############

plot-pf: 
	@./scripts/plot_pf_dist.py \
		benchmark-one/results-prefilter/psc/s12.0.txt \
		benchmark-one/results-prefilter/psc/s10.0.txt \
		benchmark-one/results-prefilter/psc/s7.5.txt
