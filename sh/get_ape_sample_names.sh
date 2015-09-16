#!/bin/bash

DIR=/home/gm114/data/great_ape_genome_project

for species in `cat $DIR/species.txt`;
do
    echo $species
    gunzip -c $DIR/VCFs/SNPs/$species.chr6_mhc.vcf.gz | grep CHR | perl -n -e '@a = split; print join "\n", @a[9..$#a]' > $DIR/samples.$species.txt
done

