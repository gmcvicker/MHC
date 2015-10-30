#!/bin/bash

# SPECIES=`sed -n ${LSB_JOBINDEX}p $SPECIES_FILE`

DIR=/home/gm114/data/great_ape_genome_project
SPECIES_FILE=$DIR/species.txt

cat $SPECIES_FILE | while read SPECIES
do
    echo $SPECIES >&2
    VCF=$DIR/VCFs/SNPs/$SPECIES.vcf.gz

    CHR6_VCF=$DIR/VCFs/SNPs/$SPECIES.chr6_mhc.vcf.gz


    gunzip -c $VCF | awk '($1 ~ /^#/ || ($1 == "chr6" && $2 > 28000000 && $2 < 34000000)) {print $0}' | gzip > $CHR6_VCF
    
    # tabix -p vcf -f $CHR6_VCF

    # if [ "$?" -ne "0" ]; then
    # 	echo "tabix failed" >&2
    # 	exit 1
    # fi
done


