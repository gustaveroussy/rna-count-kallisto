#!/bin/bash

set -e

python3.7 ../scripts/prepare_config.py genome/transcriptome.fasta --kallisto_index_extra "--kmer-size=5" --kallisto_quant_extra ""

python3.7 ../scripts/prepare_design.py reads/

snakemake -s ../Snakefile --use-conda
snakemake -s ../Snakefile --use-conda --report
