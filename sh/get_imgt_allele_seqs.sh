#!/bin/bash

GENE_FILE=$HOME/data/IMGT/gene_names.txt
SCRIPT=$HOME/proj/MHC/python/get_imgt_alleles.py

cat $GENE_FILE | while read gene;
do
    OUT_DIR=$HOME/data/IMGT/extracted/$gene
    mkdir -p $OUT_DIR
    echo $gene >&2
    python $SCRIPT --gene $gene $OUT_DIR    
done
