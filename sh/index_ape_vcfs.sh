#!/bin/bash

# SPECIES=`sed -n ${LSB_JOBINDEX}p $SPECIES_FILE`

DIR=/home/gm114/data/great_ape_genome_project
SPECIES_FILE=$DIR/species.txt

cat $SPECIES_FILE | while read SPECIES
do
    echo $SPECIES >&2
    VCF=$DIR/VCFs/SNPs/$SPECIES.vcf.gz

    tabix -p vcf -f $VCF

    if [ "$?" -ne "0" ]; then
	echo "tabix failed" >&2
	exit 1
    fi
done


