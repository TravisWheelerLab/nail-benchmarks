MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
.DEFAULT_GOAL := none
# a recipe that dies partway leaves no target behind. it reaches the
# decompression step's output and pfam.hmm, not the downloads: those land in a
# sidecar (pfam.sto.gz, swissprot.tgz) or under a phony rule (mgnify), and make
# only deletes the target named by the failing rule
.DELETE_ON_ERROR:

# detect platform + best x86_64 simd
OS   := $(shell uname -s)
ARCH := $(shell uname -m)

ifeq ($(OS),Darwin)
  PLATFORM := macos-universal
else ifeq ($(OS),Linux)
  ifeq ($(ARCH),aarch64)
    PLATFORM := linux-arm64
  else ifeq ($(ARCH),x86_64)
    CPUFLAGS := $(shell { lscpu 2>/dev/null || cat /proc/cpuinfo 2>/dev/null; } \
                        | awk -F: 'tolower($$1) ~ /flags/ { print tolower($$2); exit }')
    # only the three mmseqs2 ships a build for: avx without avx2 implies sse4.2,
    # and sse2 is part of the x86_64 base ISA, so a flag list we cannot read
    # costs the faster build rather than the build
    ifneq (,$(findstring avx2,$(CPUFLAGS)))
      SIMD := avx2
    else ifneq (,$(findstring sse4_1,$(CPUFLAGS)))
      SIMD := sse4.1
    else
      SIMD := sse2
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
PFAM_HMM := $(DATA_DIR)/pfam.hmm
HMMBUILD_CPU ?= 4

SWISSPROT_URL := https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2023_05/knowledgebase/uniprot_sprot-only2023_05.tar.gz
SWISSPROT_TGZ := $(DATA_DIR)/swissprot.tgz
SWISSPROT_DIR := $(DATA_DIR)/uniprot_sprot/
SWISSPROT_FA_GZ := $(SWISSPROT_DIR)/uniprot_sprot.fasta.gz 
SWISSPROT_FA := $(DATA_DIR)/swissprot.fa

MGY_URL = https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2024_04
MGY_DIR = $(DATA_DIR)/mgnify
MGY_SHARDS = 25

.PHONY: none
none:
	@true

$(DATA_DIR):
	@mkdir -p $@

$(MGY_DIR): | $(DATA_DIR)
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

# derived rather than downloaded, and hmmbuild has to exist first, which is why
# `data` does not build it
$(PFAM_HMM): $(PFAM_STO)
	$(call need_tool,$(HMMBUILD),hmmer)
	$(HMMBUILD) --cpu $(HMMBUILD_CPU) $@ $(PFAM_STO)

.PHONY: pfam pfam-hmm swissprot mgnify
pfam: $(PFAM_STO)

pfam-hmm: $(PFAM_HMM)

swissprot: $(SWISSPROT_FA)

mgnify: | $(MGY_DIR)
	@set -e; \
	if command -v aria2c >/dev/null 2>&1; then \
	  for i in $$(seq 1 $(MGY_SHARDS)); do \
	    aria2c -x8 -s8 -d $(MGY_DIR) $(MGY_URL)/mgy_proteins_$${i}.fa.gz; \
	  done; \
	else \
	  for i in $$(seq 1 $(MGY_SHARDS)); do \
	    wget -O $(MGY_DIR)/mgy_proteins_$${i}.fa.gz $(MGY_URL)/mgy_proteins_$${i}.fa.gz; \
	  done; \
	fi

.PHONY: data
data: pfam swissprot mgnify

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

# stop at the platform rather than handing wget the string "none"
need_url = @test "$(2)" != none || { echo "no $(1) binary release for $(PLATFORM)" >&2; exit 1; }

# $(1) is the binary, $(2) the target that installs it
need_tool = @test -x $(1) || { echo "no $(notdir $(1)); run make $(2)" >&2; exit 1; }

TOOL_DIR := $(MAKEFILE_DIR)/tools
TOOL_BIN := $(TOOL_DIR)/bin

$(TOOL_BIN):
	@mkdir -p $@

NAIL        := $(TOOL_BIN)/nail
PHMMER      := $(TOOL_BIN)/phmmer
HMMSEARCH   := $(TOOL_BIN)/hmmsearch
ESL_SEQSTAT := $(TOOL_BIN)/esl-seqstat
ESL_ALISTAT := $(TOOL_BIN)/esl-alistat
HMMSTAT     := $(TOOL_BIN)/hmmstat
PROFMARK    := $(TOOL_BIN)/create-profmark
HMMBUILD    := $(TOOL_BIN)/hmmbuild
HMMEMIT     := $(TOOL_BIN)/hmmemit
MMSEQS      := $(TOOL_BIN)/mmseqs
LASTAL      := $(TOOL_BIN)/lastal
LASTDB      := $(TOOL_BIN)/lastdb
BLASTP      := $(TOOL_BIN)/blastp
PSIBLAST    := $(TOOL_BIN)/psiblast
MAKEBLASTDB := $(TOOL_BIN)/makeblastdb
DIAMOND     := $(TOOL_BIN)/diamond

.PHONY: nail hmmer mmseqs blast last diamond

NAIL_SRC_URL  := https://github.com/TravisWheelerLab/nail/archive/refs/tags/nail-v0.7.1.tar.gz
NAIL_SRC_TGZ  := $(TOOL_DIR)/nail.tgz
NAIL_SRC_DIR  := $(TOOL_DIR)/nail
NAIL_TGT_DIR  := $(NAIL_SRC_DIR)/target
# the build passes --target-dir explicitly: CARGO_TARGET_DIR in the environment
# would otherwise move the binary and leave the symlink dangling, since ln -s
# does not check
nail: $(TOOL_BIN)
	@wget -O $(NAIL_SRC_TGZ) $(NAIL_SRC_URL)
	@mkdir -p $(NAIL_SRC_DIR)
	@tar --strip-components=1 -xzf $(NAIL_SRC_TGZ) -C $(NAIL_SRC_DIR)
	@cd $(NAIL_SRC_DIR) && cargo build --release -p nail --target-dir $(NAIL_TGT_DIR)
	@rm $(NAIL_SRC_TGZ)
	@ln -sf $(NAIL_TGT_DIR)/release/nail $(NAIL)

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
	@ln -sf $(HMMER_BIN_DIR)/esl-alistat $(ESL_ALISTAT)
	@ln -sf $(HMMER_BIN_DIR)/hmmstat $(HMMSTAT)
	@ln -sf $(HMMER_BIN_DIR)/hmmbuild $(HMMBUILD)
	@ln -sf $(HMMER_BIN_DIR)/hmmemit $(HMMEMIT)
	@ln -sf $(HMMER_SRC_DIR)/profmark/create-profmark $(PROFMARK)
	@rm $(HMMER_SRC_TGZ)

MMSEQS_BIN_TGZ := $(TOOL_DIR)/mmseqs.tgz
MMSEQS_DIR     := $(TOOL_DIR)/mmseqs
MMSEQS_BIN     := $(MMSEQS_DIR)/bin/mmseqs
mmseqs: $(TOOL_BIN)
	$(call need_url,mmseqs,$(MMSEQS_BIN_URL))
	@wget -O $(MMSEQS_BIN_TGZ) $(MMSEQS_BIN_URL)
	@mkdir -p $(MMSEQS_DIR)
	@tar --strip-components=1 -xzf $(MMSEQS_BIN_TGZ) -C $(MMSEQS_DIR)
	@rm $(MMSEQS_BIN_TGZ)
	@ln -sf $(MMSEQS_BIN) $(MMSEQS)

LAST_SRC_URL := https://gitlab.com/mcfrith/last/-/archive/1642/last-1642.tar.gz
LAST_SRC_TGZ := $(TOOL_DIR)/last.tgz
LAST_SRC_DIR := $(TOOL_DIR)/last
LAST_BIN_DIR := $(LAST_SRC_DIR)/bin
last: $(TOOL_BIN)
	@wget -O $(LAST_SRC_TGZ) $(LAST_SRC_URL)
	@mkdir -p $(LAST_SRC_DIR)
	@tar --strip-components=1 -xzf $(LAST_SRC_TGZ) -C $(LAST_SRC_DIR)
	@cd $(LAST_SRC_DIR) && make
	@rm $(LAST_SRC_TGZ)
	@ln -sf $(LAST_BIN_DIR)/lastal $(LASTAL)
	@ln -sf $(LAST_BIN_DIR)/lastdb $(LASTDB)

BLAST_BIN_TGZ := $(TOOL_DIR)/blast.tgz
BLAST_DIR     := $(TOOL_DIR)/blast
BLAST_BIN_DIR := $(BLAST_DIR)/bin
blast: $(TOOL_BIN)
	$(call need_url,blast,$(BLAST_BIN_URL))
	@wget -O $(BLAST_BIN_TGZ) $(BLAST_BIN_URL)
	@mkdir -p $(BLAST_DIR)
	@tar --strip-components=1 -xzf $(BLAST_BIN_TGZ) -C $(BLAST_DIR)
	@rm $(BLAST_BIN_TGZ)
	@ln -sf $(BLAST_BIN_DIR)/blastp $(BLASTP)
	@ln -sf $(BLAST_BIN_DIR)/psiblast $(PSIBLAST)
	@ln -sf $(BLAST_BIN_DIR)/makeblastdb $(MAKEBLASTDB)

DIAMOND_BIN_TGZ := $(TOOL_DIR)/diamond.tgz
diamond: $(TOOL_BIN)
	$(call need_url,diamond,$(DIAMOND_BIN_URL))
	@wget -O $(DIAMOND_BIN_TGZ) $(DIAMOND_BIN_URL)
	@tar -xzf $(DIAMOND_BIN_TGZ) -C $(TOOL_BIN)
	@rm $(DIAMOND_BIN_TGZ)

# built from source, so every platform can have them
TOOLS := nail hmmer last
# the rest ship binaries, and not for every platform
ifneq ($(MMSEQS_BIN_URL),none)
  TOOLS += mmseqs
endif
ifneq ($(BLAST_BIN_URL),none)
  TOOLS += blast
endif
ifneq ($(DIAMOND_BIN_URL),none)
  TOOLS += diamond
endif

.PHONY: tools
tools: $(TOOLS)

# name:help-flag, matching what util::tools runs before handing back a path
CHECK_TOOLS := nail:-h hmmsearch:-h phmmer:-h hmmbuild:-h hmmemit:-h \
               esl-seqstat:-h create-profmark:-h mmseqs:-h \
               blastp:-h psiblast:-h makeblastdb:-h \
               lastal:-h lastdb:-h diamond:--help

CHECK_DATA := pfam.sto pfam.hmm swissprot.fa mgnify mgy-cutoffs.tbl long-seqs

.PHONY: check
check:
	@fail=0; \
	echo "tools  $(TOOL_BIN)"; \
	for spec in $(CHECK_TOOLS); do \
	  name=$${spec%%:*}; flag=$${spec#*:}; bin=$(TOOL_BIN)/$$name; note=; \
	  if [ -L "$$bin" ] && [ ! -e "$$bin" ]; then \
	    note="dangling symlink -> $$(readlink $$bin)"; \
	  elif [ ! -e "$$bin" ]; then \
	    note="not installed"; \
	  elif ! "$$bin" $$flag >/dev/null 2>&1; then \
	    note="$$flag exited nonzero"; \
	  fi; \
	  if [ -n "$$note" ]; then \
	    fail=1; printf '  ✘ %-18s %s\n' "$$name" "$$note"; \
	  else printf '  ✔ %s\n' "$$name"; fi; \
	done; \
	echo; \
	echo "data  $(DATA_DIR)"; \
	for item in $(CHECK_DATA); do \
	  path=$(DATA_DIR)/$$item; note=; \
	  if [ ! -e "$$path" ]; then \
	    note="missing"; \
	    if [ "$$item" = pfam.hmm ]; then \
	      note="missing; run make pfam-hmm"; \
	    fi; \
	  elif [ "$$item" = mgnify ]; then \
	    n=$$(ls $(MGY_DIR)/*.fa $(MGY_DIR)/*.fasta 2>/dev/null | wc -l | tr -d ' '); \
	    if [ "$$n" = 0 ]; then note="no .fa/.fasta in it; the downloads are .gz"; fi; \
	  fi; \
	  if [ -n "$$note" ]; then \
	    fail=1; printf '  ✘ %-18s %s\n' "$$item" "$$note"; \
	  else printf '  ✔ %s\n' "$$item"; fi; \
	done; \
	echo; \
	if [ $$fail = 0 ]; then \
	  echo "all present"; \
	else \
	  echo "something is missing: make tools, make data"; \
	fi; \
	exit $$fail

# check says a file is there; validate looks inside. each file goes to the tool
# that reads its format, and the count that tool reports is printed. mgnify is
# the exception: it is the one input that runs to hundreds of gigabytes, and a
# fasta cannot be validated by reading anyway -- any prefix of one is a valid
# fasta -- so it is measured instead. mgy-cutoffs.tbl has no reader, so its rows
# are counted directly
#
# nothing here is compared against anything: no target knows how many records a
# file holds. the numbers are printed so a wrong one is visible
.PHONY: validate
validate:
	$(call need_tool,$(ESL_ALISTAT),hmmer)
	$(call need_tool,$(ESL_SEQSTAT),hmmer)
	$(call need_tool,$(HMMSTAT),hmmer)
	@fail=0; seen=0; \
	ok() { printf '  ✔ %-20s %s\n' "$$1" "$$2"; seen=$$((seen + 1)); }; \
	bad() { printf '  ✘ %-20s %s\n' "$$1" "$$2"; seen=$$((seen + 1)); fail=$$((fail + 1)); }; \
	records() { printf '%s\n' "$$1" | awk '!/^#/ && NF { n++ } END { print n+0 }'; }; \
	nseq() { printf '%s\n' "$$1" | awk '/^Number of sequences:/ { print $$4 }'; }; \
	echo "data  $(DATA_DIR)"; \
	if [ ! -e $(PFAM_STO) ]; then bad pfam.sto missing; \
	elif out=$$($(ESL_ALISTAT) -1 --informat stockholm $(PFAM_STO) 2>&1); then \
	  ok pfam.sto "$$(records "$$out") alignments"; \
	else bad pfam.sto "esl-alistat did not finish"; fi; \
	if [ ! -e $(PFAM_HMM) ]; then bad pfam.hmm "missing; run make pfam-hmm"; \
	elif out=$$($(HMMSTAT) $(PFAM_HMM) 2>&1); then \
	  ok pfam.hmm "$$(records "$$out") models"; \
	else bad pfam.hmm "hmmstat did not finish"; fi; \
	if [ ! -e $(SWISSPROT_FA) ]; then bad swissprot.fa missing; \
	elif out=$$($(ESL_SEQSTAT) $(SWISSPROT_FA) 2>&1); then \
	  ok swissprot.fa "$$(nseq "$$out") sequences"; \
	else bad swissprot.fa "esl-seqstat did not finish"; fi; \
	if [ ! -d $(MGY_DIR) ]; then bad mgnify missing; \
	else \
	  n=$$(ls $(MGY_DIR) | wc -l | tr -d ' '); \
	  if [ "$$n" = 0 ]; then bad mgnify empty; \
	  else ok mgnify "$$(du -sh $(MGY_DIR) | cut -f1) in $$n files"; fi; \
	fi; \
	if [ ! -d $(DATA_DIR)/long-seqs ]; then bad long-seqs missing; \
	else \
	  n=0; total=0; broke=; \
	  for f in $(DATA_DIR)/long-seqs/*/*.fa; do \
	    [ -e "$$f" ] || continue; \
	    n=$$((n + 1)); \
	    if out=$$($(ESL_SEQSTAT) "$$f" 2>&1); then total=$$((total + $$(nseq "$$out"))); \
	    else broke=$$(basename "$$f"); fi; \
	  done; \
	  if [ -n "$$broke" ]; then bad long-seqs "esl-seqstat did not finish on $$broke"; \
	  elif [ $$n = 0 ]; then bad long-seqs "no .fa in it"; \
	  else ok long-seqs "$$total sequences in $$n files"; fi; \
	fi; \
	if [ ! -e $(DATA_DIR)/mgy-cutoffs.tbl ]; then bad mgy-cutoffs.tbl missing; \
	else ok mgy-cutoffs.tbl \
	  "$$(awk '!/^#/ && NF { n++ } END { print n+0 }' $(DATA_DIR)/mgy-cutoffs.tbl) rows"; fi; \
	echo; \
	if [ $$fail = 0 ]; then echo "$$seen checks, no problems"; \
	else echo "$$seen checks, $$fail with problems"; exit 1; fi

.PHONY: clean
clean:
	rm -rf $(TOOL_DIR)
