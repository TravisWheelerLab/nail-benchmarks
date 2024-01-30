# nail-benchmarks
A collection of reproducible benchmarks for [nail](https://github.com/travisWheelerLab/nail)

## Requirements

You'll need the following available on your system path:

- [nail](https://github.com/travisWheelerLab/nail)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [HMMER3](http://hmmer.org/)
- easel (comes with HMMER3 distributions)
- the create-profmark binary (comes with HMMER3 distributions)

## Download sequence data
This benchmark was originally run using Pfam version `36.0` and Swissprot `release-2023_05`

To download the data, you can run
    
    $ ./scripts/download-data.sh

which will place Pfam seed alignments & Swissprot sequences in the `data/` directory:

```
$ tree data
data
├── long-seq
│   ├── query
│   │   ├── 1.query.fa
│   │   ├── 2.query.fa
│   │   ├── 3.query.fa
│   │   ├── 4.query.fa
│   │   ├── 5.query.fa
│   │   └── 6.query.fa
│   └── target
│       ├── 1.target.fa
│       ├── 2.target.fa
│       ├── 3.target.fa
│       ├── 4.target.fa
│       ├── 5.target.fa
│       └── 6.target.fa
├── pfam.sto
├── uniprot.tar.gz
├── uniprot_sprot.dat.gz
├── uniprot_sprot.fasta
├── uniprot_sprot.fasta.ssi
├── uniprot_sprot.xml.gz
└── uniprot_sprot_varsplic.fasta.gz
```

## Build the benchmark

To build the benchmark, run

    $ ./scripts/build-benchmark.sh

## Run the benchmark

To run the benchmark, run

    $ ./scripts/run-all.sh

## Produce plots

To produce the plots, run

    $ python ./scripts/plots.py ./benchmark/
