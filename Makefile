MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
.DEFAULT_GOAL := none

DATA_DIR := $(MAKEFILE_DIR)/data

PFAM_URL := https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam36.0/Pfam-A.seed.gz
PFAM_GZ := $(DATA_DIR)/pfam.sto.gz
PFAM_STO:= $(DATA_DIR)/pfam.sto

SWISSPROT_URL := https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2023_05/knowledgebase/uniprot_sprot-only2023_05.tar.gz
SWISSPROT_TGZ := $(DATA_DIR)/swissprot.tgz
SWISSPROT_DIR := $(DATA_DIR)/uniprot_sprot/
SWISSPROT_FA_GZ := $(SWISSPROT_DIR)/uniprot_sprot.fasta.gz 
SWISSPROT_FA := $(DATA_DIR)/swissprot.fasta

.PHONY: none
none:
	@true

$(DATA_DIR):
	@mkdir -p $@

$(SWISSPROT_FA): | $(DATA_DIR)
	@wget -O $(SWISSPROT_TGZ) $(SWISSPROT_URL)
	@mkdir -p $(SWISSPROT_DIR)
	@tar -xzf $(SWISSPROT_TGZ) -C $(SWISSPROT_DIR)
	@gunzip -c $(SWISSPROT_FA_GZ) > $(SWISSPROT_FA)
	@rm -rf $(SWISSPROT_TGZ)
	@rm -rf $(SWISSPROT_DIR)

$(PFAM_STO): | $(DATA_DIR)
	@wget -O $(PFAM_GZ) $(PFAM_URL)
	@gunzip $(PFAM_GZ)

DATA := $(PFAM_STO) $(SWISSPROT_FA)
setup: $(DATA)

####################################
####################################
####################################

TOOLS_DIR := $(MAKEFILE_DIR)/tools
TOOL_BIN := $(TOOLS_DIR)/bin

$(TOOL_BIN):
	@mkdir -p $@

NAIL        := $(TOOL_BIN)/nail
PHMMER      := $(TOOL_BIN)/phmmer
HMMSEARCH   := $(TOOL_BIN)/hmmsearch
ESL_SEQSTAT := $(TOOL_BIN)/esl-seqstat
PROFMARK    := $(TOOL_BIN)/create-profmark
MMSEQS      := $(TOOL_BIN)/mmseqs
BLAST       := $(TOOL_BIN)/blastp
LAST        := $(TOOL_BIN)/lastal

.PHONY: nail hmmer mmseqs blast last

NAIL_SRC_URL  := https://github.com/TravisWheelerLab/nail/archive/refs/tags/nail-v0.4.0.tar.gz NAIL_SRC_TGZ  := $(TOOLS_DIR)/nail.tgz
nail: $(TOOL_BIN)
	@echo TODO: retrieve/build nail
	@False	

HMMER_SRC_URL := http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz
HMMER_SRC_TGZ := $(TOOLS_DIR)/hmmer.tgz
HMMER_SRC_DIR := $(TOOLS_DIR)/hmmer
HMMER_BIN_DIR := $(HMMER_SRC)/bin/
hmmer: $(TOOL_BIN)
	@wget -O $(HMMER_SRC_TGZ) $(HMMER_SRC_URL)
	@mkdir -p $(HMMER_SRC_DIR)
	@tar --strip-components=1 -xzf $(HMMER_SRC_TGZ) -C $(HMMER_SRC_DIR)
	@cd $(HMMER_SRC_DIR) && \
		./configure && \
		make install prefix=$(HMMER_SRC_DIR) && \
		cd easel && \
		make install prefix=$(HMMER_SRC_DIR)
	@ln -sf $(HMMER_BIN_DIR)/hmmsearch $(HMMSEARCH)
	@ln -sf $(HMMER_BIN_DIR)/phmmer $(PHMMER)
	@ln -sf $(HMMER_BIN_DIR)/esl-seqstat $(ESL_SEQSTAT)
	@ln -sf $(HMMER_SRC_DIR)/profmark/create-profmark $(PROFMARK)
	@rm $(HMMER_SRC_TGZ)

MMSEQS_BIN_URL := https://github.com/soedinglab/MMseqs2/releases/download/18-8cc5c/mmseqs-osx-universal.tar.gz
MMSEQS_BIN_TGZ := $(TOOLS_DIR)/mmseqs.tgz
MMSEQS_BIN_DIR := $(TOOLS_DIR)/mmseqs
MMSEQS_BIN     := $(MMSEQS_BIN_DIR)/bin/mmseqs
mmseqs: $(TOOL_BIN)
	@wget -O $(MMSEQS_BIN_TGZ) $(MMSEQS_BIN_URL)
	@mkdir -p $(MMSEQS_BIN_DIR)
	@tar --strip-components=1 -xzf $(MMSEQS_BIN_TGZ) -C $(MMSEQS_BIN_DIR)
	@rm $(MMSEQS_BIN_TGZ)
	@ln -sf $(MMSEQS_BIN) $(MMSEQS)

# BLAST_BIN_URL := https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-universal-macosx.tar.gz
# BLAST_BIN_TGZ := $(TOOLS_DIR)/blast.tgz
# BLAST_BIN_DIR := $(TOOLS_DIR)/blast
# BLAST_BIN_DIR := $(BLAST_BIN_DIR)/bin
# blast: $(TOOL_BIN)
# 	@wget -O $(BLAST_BIN_TGZ) $(BLAST_BIN_URL)
# 	@mkdir -p $(BLAST_BIN_DIR)
# 	@tar --strip-components=1 -xzf $(BLAST_BIN_TGZ) -C $(BLAST_BIN_DIR)
# 	@rm $(BLAST_BIN_TGZ)
# 	@ln -sf $(BLAST_BIN) $(BLAST)

# LAST_SRC_URL := https://gitlab.com/mcfrith/last/-/archive/1642/last-1642.tar.gz
# LAST_SRC_TGZ := $(TOOLS_DIR)/last.tgz
# LAST_SRC_DIR := $(TOOLS_DIR)/last
# LAST_BIN_DIR := $(LAST_SRC)/bin
# last: $(TOOL_BIN)
# 	@wget -O $(LAST_SRC_TGZ) $(LAST_SRC_URL)
# 	@mkdir -p $(LAST_SRC)
# 	@tar --strip-components=1 -xzf $(LAST_SRC_TGZ) -C $(LAST_SRC)
# 	@cd $(LAST_SRC) && make
# 	@rm $(LAST_SRC_TGZ)
# 	@ln -sf $(LAST_BIN) $(LAST)
