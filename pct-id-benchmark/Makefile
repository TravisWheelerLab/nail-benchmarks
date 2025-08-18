PM_BIN := create-profmark
SPLIT_BIN := ./scripts/split-decoys

PM_DIR := profmark
PM_NAME := $(PM_DIR)/benchmark

SRC_MSA := $(realpath ../data/pfam.sto)
# SRC_MSA := $(realpath ../data/pfam.100.sto)
SRC_FA := $(realpath ../data/uniprot_sprot.fasta)
SRC_FA_IDX := $(FA).ssi

TRAIN_MSA := $(PM_NAME).train.msa
TEST_FA := $(PM_NAME).test.fa

POS_FA := $(PM_DIR)/positives.fa
DEC_FA := $(PM_DIR)/decoys.fa
SRC_STO := $(PM_DIR)/source.sto
QUERY_STO := $(PM_DIR)/query.sto
QUERY_HMM := $(PM_DIR)/query.hmm
QUERY_CONS := $(PM_DIR)/query.cons.fa

PM_OUT := $(POS_FA) $(DEC_FA) $(QUERY_STO) $(QUERY_HMM) $(QUERY_CONS)
PM_REMOVE := $(TEST_FA) $(PM_NAME).neg $(PM_NAME).pos $(PM_NAME).tbl

SEED := 67779
N_DECOYS := 5000000
# N_DECOYS := 1000
TRAIN_TEST_ID := 0.5
MIN_TEST := 10
MAX_TEST := 30

BUILD_BIN := ./scripts/pid
BUILD_ARGS := $(POS_FA) $(DEC_FA) $(QUERY_STO) $(SRC_STO)
DECOY_RATIO := 100

.PHONY: 100 1000 benchmark

$(FA_IDX): 
	@echo indexing \"$(FA)\"...
	@esl-sfetch --index $(FA) >/dev/null

$(PM_OUT): $(FA_IDX)
	@echo running profmark...
	@mkdir -p $(PM_DIR)
	@time $(PM_BIN) \
		-S $(SEED) \
		-N $(N_DECOYS) \
		-1 $(TRAIN_TEST_ID) \
		--cluster \
		--mintest $(MIN_TEST) \
		--maxtest $(MAX_TEST) \
		$(PM_NAME) $(SRC_MSA) $(SRC_FA)
	
	@$(SPLIT_BIN) $(TEST_FA)
	@mv $(TRAIN_MSA) $(QUERY_STO)
	@hmmbuild --cpu 8 $(QUERY_HMM) $(QUERY_STO) > /dev/null
	@hmmemit -c $(QUERY_HMM) > $(QUERY_CONS)
	@ln -sf $(SRC_MSA) $(SRC_STO)
	@rm $(PM_REMOVE)
		
profmark: $(PM_OUT)

10: profmark
	@mkdir -p b10
	@ln -sf ../$(QUERY_HMM)  b10/query.hmm
	@ln -sf ../$(QUERY_CONS) b10/query.cons.fa
	@time $(BUILD_BIN) b10/ $(BUILD_ARGS) 10 $(DECOY_RATIO)

100: profmark
	@mkdir -p b100
	@ln -sf ../$(QUERY_HMM)  b100/query.hmm
	@ln -sf ../$(QUERY_CONS) b100/query.cons.fa
	@time $(BUILD_BIN) b100/ $(BUILD_ARGS) 100 $(DECOY_RATIO)

benchmark: profmark
	@mkdir -p benchmark
	@ln -sf ../$(QUERY_HMM)  benchmark/query.hmm
	@ln -sf ../$(QUERY_CONS) benchmark/query.cons.fa
	@time $(BUILD_BIN) benchmark/ $(BUILD_ARGS) 1000 $(DECOY_RATIO)

clean:
	@rm -rf $(PM_DIR)/ b10/ b100/ benchmark/
