#!/bin/bash
#BSUB -n 1
#BSUB -W 0:30
#BSUB -J jobArray[1-261]
#BSUB -o /home/gm114/lsf/%J.%I.out
#BSUB -e /home/gm114/lsf/%J.%I.err
#BSUB -q short

#
# Re-map SGDP contigs that mapped to any of
# the MHC haplotypes, but this time map them
# all to the reference genome so in same coord system...
#

REF_GENOME=/home/gm114/data/GRCh38/ref_chroms/all_ref_chroms.fa.gz
CONTIG_DIR=/home/gm114/data/SGDP/mag_mhc_mapped

SAMPLES_FILE=/home/gm114/data/SGDP/C-team-mag/261.samples

SAMPLE_ID=`sed -n ${LSB_JOBINDEX}p $SAMPLES_FILE | awk '{print $2}'`
CONTIG_FILE=$CONTIG_DIR/$SAMPLE_ID.fq.gz

OUTPUT_DIR=/home/gm114/data/SGDP/mag_mhc_mapped/
mkdir -p $OUTPUT_DIR

bwa mem $REF_GENOME $CONTIG_FILE | gzip > $OUTPUT_DIR/$SAMPLE_ID.sam.gz

if [ "$?" -ne "0" ]; then
    echo "bwa mem failed";
    exit 1
fi

