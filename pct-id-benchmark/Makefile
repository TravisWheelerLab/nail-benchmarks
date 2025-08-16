BUILD := ./scripts/build-benchmark.sh
FASTA := ../data/uniprot_sprot.fasta

.PHONY: 100
100:
	@rm -rf b100/
	@$(BUILD) ../data/pfam.100.sto $(FASTA) 0 b100/

.PHONY: 1000
1000:
	@rm -rf b1000/
	@$(BUILD) ../data/pfam.1000.sto $(FASTA) 0 b1000/

.PHONY: benchmark
benchmark:
	@rm -rf benchmark/
	@$(BUILD) ../data/pfam.sto $(FASTA) 0 benchmark/
