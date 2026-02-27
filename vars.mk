ROOT  := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
TOOLS := $(ROOT)/tools
UTIL  := $(ROOT)/util
DATA  := $(ROOT)/data

RUST     := $(UTIL)/rust
RUST_BIN := $(RUST)/target/release

BIOIO_SRC := $(wildcard $(RUST)/bioio/src/*.rs)

# find all the <bin>.rs files
BINS := $(patsubst %.rs,%,$(notdir $(wildcard $(RUST)/bm/src/bin/*.rs)))

# create uppercase make vars for each binary
$(foreach b,$(BINS),$(eval $(shell echo $(b) | tr a-z- A-Z_) := $(RUST_BIN)/$(b)))

# create a target for every binary
$(BINS): %: $(RUST_BIN)/% $(BIOIO_SRC)
	@:

# actual binary builds
$(RUST_BIN)/%: $(RUST)/bm/src/bin/%.rs
	@cd $(RUST) && cargo build --release --bin $*

RUN_BINS := $(addprefix bin-,$(BINS))

.PHONY: $(RUN_BINS)

$(RUN_BINS):
	@$(RUST_BIN)/$(@:bin-%=%)

FASTA_BALANCE := $(UTIL)/scripts/fastabalance.py
FASTA_BIN     := $(UTIL)/scripts/fastabin.py
FASTA_SAMPLE  := $(UTIL)/scripts/fastasample.py

HMM_BALANCE   := $(UTIL)/scripts/hmmbalance.py
HMM_BIN       := $(UTIL)/scripts/hmmbin.py

STO_STAMPLE   := $(UTIL)/scripts/sample-sto.sh

# hmmer util
PROFMARK   := $(TOOLS)/bin/create-profmark
HMMBUILD   := $(TOOLS)/bin/hmmbuild
HMMEMIT    := $(TOOLS)/bin/hmmemit



