MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
.DEFAULT_GOAL := none

# detect platform + best x86_64 simd
OS   := $(shell uname -s)
ARCH := $(shell uname -m)

ifeq ($(OS),Darwin)
  PLATFORM := macos-universal
else ifeq ($(OS),Linux)
  ifeq ($(ARCH),aarch64)
    PLATFORM := linux-arm64
  else ifeq ($(ARCH),x86_64)
    CPUFLAGS := $(shell lscpu 2>/dev/null | awk -F: '/Flags|flags/{print $$2}' | tr A-Z a-z); \
                if [ -z "$$CPUFLAGS" ]; then CPUFLAGS=$$(grep -im1 '^flags' /proc/cpuinfo | cut -d: -f2); fi; \
                echo $$CPUFLAGS
    ifneq (,$(findstring avx2,$(CPUFLAGS)))
      SIMD := avx2
    else ifneq (,$(findstring avx,$(CPUFLAGS)))
      SIMD := avx
    else ifneq (,$(findstring sse4_1,$(CPUFLAGS)))
      SIMD := sse4.1
    else ifneq (,$(findstring sse2,$(CPUFLAGS)))
      SIMD := sse2
    else
      SIMD := baseline
    endif
    PLATFORM := linux-x86_64-$(SIMD)
  else
    PLATFORM := linux-$(ARCH)
  endif
else
  PLATFORM := unknown
endif

print-platform:
	@echo $(PLATFORM)

DATA_DIR := $(MAKEFILE_DIR)/data

PFAM_URL := https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam36.0/Pfam-A.seed.gz
PFAM_GZ := $(DATA_DIR)/pfam.sto.gz
PFAM_STO:= $(DATA_DIR)/pfam.sto

SWISSPROT_URL := https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2023_05/knowledgebase/uniprot_sprot-only2023_05.tar.gz
SWISSPROT_TGZ := $(DATA_DIR)/swissprot.tgz
SWISSPROT_DIR := $(DATA_DIR)/uniprot_sprot/
SWISSPROT_FA_GZ := $(SWISSPROT_DIR)/uniprot_sprot.fasta.gz 
SWISSPROT_FA := $(DATA_DIR)/swissprot.fa

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

ifeq ($(PLATFORM),linux-arm64)
  MMSEQS_BIN_URL  := https://github.com/soedinglab/MMseqs2/releases/download/18-8cc5c/mmseqs-linux-arm64.tar.gz
  BLAST_BIN_URL   := https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-aarch64-linux.tar.gz
  DIAMOND_BIN_URL := none
else ifeq ($(PLATFORM),linux-x86_64-avx2)
  MMSEQS_BIN_URL  := https://github.com/soedinglab/MMseqs2/releases/download/18-8cc5c/mmseqs-linux-avx2.tar.gz
  BLAST_BIN_URL   := https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
  DIAMOND_BIN_URL := https://github.com/bbuchfink/diamond/releases/download/v2.1.13/diamond-linux64.tar.gz
else ifeq ($(PLATFORM),linux-x86_64-sse2)
  MMSEQS_BIN_URL  := https://github.com/soedinglab/MMseqs2/releases/download/18-8cc5c/mmseqs-linux-sse2.tar.gz
  BLAST_BIN_URL   := https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
  DIAMOND_BIN_URL := https://github.com/bbuchfink/diamond/releases/download/v2.1.13/diamond-linux64.tar.gz
else ifeq ($(PLATFORM),linux-x86_64-sse4.1)
  MMSEQS_BIN_URL  := https://github.com/soedinglab/MMseqs2/releases/download/18-8cc5c/mmseqs-linux-sse41.tar.gz
  BLAST_BIN_URL   := https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
  DIAMOND_BIN_URL := https://github.com/bbuchfink/diamond/releases/download/v2.1.13/diamond-linux64.tar.gz
else ifeq ($(PLATFORM),macos-universal)
  MMSEQS_BIN_URL  := https://github.com/soedinglab/MMseqs2/releases/download/18-8cc5c/mmseqs-osx-universal.tar.gz
  BLAST_BIN_URL   := https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-universal-macosx.tar.gz
  DIAMOND_BIN_URL := https://github.com/bbuchfink/diamond/releases/download/v2.1.13/diamond-macos.tar.gz
else
  MMSEQS_BIN_URL  := none
  BLAST_BIN_URL   := none
  DIAMOND_BIN_URL := none
endif

TOOL_DIR := $(MAKEFILE_DIR)/tools
TOOL_BIN := $(TOOL_DIR)/bin

$(TOOL_BIN):
	@mkdir -p $@

NAIL        := $(TOOL_BIN)/nail
PHMMER      := $(TOOL_BIN)/phmmer
HMMSEARCH   := $(TOOL_BIN)/hmmsearch
ESL_SEQSTAT := $(TOOL_BIN)/esl-seqstat
PROFMARK    := $(TOOL_BIN)/create-profmark
HMMBUILD    := $(TOOL_BIN)/hmmbuild
HMMEMIT     := $(TOOL_BIN)/hmmemit
MMSEQS      := $(TOOL_BIN)/mmseqs
BLASTP      := $(TOOL_BIN)/blastp
PSIBLAST    := $(TOOL_BIN)/psiblast
MAKEBLASTDB := $(TOOL_BIN)/makeblastdb
DIAMOND     := $(TOOL_BIN)/diamond

.PHONY: nail hmmer mmseqs blast last diamond

NAIL_SRC_URL  := https://github.com/TravisWheelerLab/nail/archive/refs/tags/nail-v0.4.0.tar.gz
NAIL_SRC_TGZ  := $(TOOL_DIR)/nail.tgz
nail: $(TOOL_BIN)
	@echo TODO: retrieve/build nail
	@False	

HMMER_SRC_URL := http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz
HMMER_SRC_TGZ := $(TOOL_DIR)/hmmer.tgz
HMMER_SRC_DIR := $(TOOL_DIR)/hmmer
HMMER_BIN_DIR := $(HMMER_SRC_DIR)/bin/
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
	@ln -sf $(HMMER_BIN_DIR)/hmmbuild $(HMMBUILD)
	@ln -sf $(HMMER_BIN_DIR)/hmmemit $(HMMEMIT)
	@ln -sf $(HMMER_SRC_DIR)/profmark/create-profmark $(PROFMARK)
	@rm $(HMMER_SRC_TGZ)

MMSEQS_BIN_TGZ := $(TOOL_DIR)/mmseqs.tgz
MMSEQS_DIR     := $(TOOL_DIR)/mmseqs
MMSEQS_BIN     := $(MMSEQS_DIR)/bin/mmseqs
mmseqs: $(TOOL_BIN)
	@wget -O $(MMSEQS_BIN_TGZ) $(MMSEQS_BIN_URL)
	@mkdir -p $(MMSEQS_DIR)
	@tar --strip-components=1 -xzf $(MMSEQS_BIN_TGZ) -C $(MMSEQS_DIR)
	@rm $(MMSEQS_BIN_TGZ)
	@ln -sf $(MMSEQS_BIN) $(MMSEQS)

BLAST_BIN_TGZ := $(TOOL_DIR)/blast.tgz
BLAST_DIR     := $(TOOL_DIR)/blast
BLAST_BIN_DIR := $(BLAST_DIR)/bin
blast: $(TOOL_BIN)
	@wget -O $(BLAST_BIN_TGZ) $(BLAST_BIN_URL)
	@mkdir -p $(BLAST_DIR)
	@tar --strip-components=1 -xzf $(BLAST_BIN_TGZ) -C $(BLAST_DIR)
	@rm $(BLAST_BIN_TGZ)
	@ln -sf $(BLAST_BIN_DIR)/blastp $(BLASTP)
	@ln -sf $(BLAST_BIN_DIR)/psiblast $(PSIBLAST)
	@ln -sf $(BLAST_BIN_DIR)/makeblastdb $(MAKEBLASTDB)

DIAMOND_BIN_TGZ := $(TOOL_DIR)/diamond.tgz
DIAMOND_BIN     := $(TOOL_DIR)/diamond
diamond: $(TOOL_BIN)
	@wget -O $(DIAMOND_BIN_TGZ) $(DIAMOND_BIN_URL)
	@tar -xzf $(DIAMOND_BIN_TGZ) -C $(TOOL_BIN)
	@rm $(DIAMOND_BIN_TGZ)

.PHONY: clean
clean:
	rm -rf $(TOOL_DIR)
