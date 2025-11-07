ROOT_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
TOOLS    := $(realpath $(ROOT_DIR)/tools)
UTIL     := $(realpath $(ROOT_DIR)/util)
DATA     := $(realpath $(ROOT_DIR)/data)

RUST     := $(UTIL)/rust
RUST_BIN := $(RUST)/target/release

BIOIO_SRC := $(wildcard $(RUST)/bioio/src/*.rs)

FASTA_BALANCE := $(UTIL)/scripts/fastabalance.py
FASTA_BIN     := $(UTIL)/scripts/fastabin.py
FASTA_SAMPLE  := $(UTIL)/scripts/fastasample.py

HMM_BALANCE   := $(UTIL)/scripts/hmmbalance.py
HMM_BIN       := $(UTIL)/scripts/hmmbin.py

STO_STAMPLE   := $(UTIL)/scripts/sample-sto.sh
