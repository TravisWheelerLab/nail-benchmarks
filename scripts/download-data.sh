#! /bin/sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

PFAM=$SCRIPT_DIR/../data/pfam.gz
UNIPROT=$SCRIPT_DIR/../data/uniprot.gz 

wget -O $PFAM https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam36.0/Pfam-A.seed.gz
#gunzip $PFAM

wget -O $UNIPROT https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2023_05/knowledgebase/uniprot_sprot-only2023_05.tar.gz
#tar -xf $UNIPROT
