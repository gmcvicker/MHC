#!/bin/bash
#BSUB -n 1
#BSUB -W 0:30
#BSUB -J jobArray[1-261]
#BSUB -o /home/gm114/lsf/%J.%I.out
#BSUB -e /home/gm114/lsf/%J.%I.err
#BSUB -q short


DATA_DIR=/home/gm114/data/SGDP/mag_mhc_mapped/

SAMPLES_FILE=/home/gm114/data/SGDP/C-team-mag/261.samples
SAMPLE_ID=`sed -n ${LSB_JOBINDEX}p $SAMPLES_FILE | awk '{print $2}'`

SAM_FILE=$DATA_DIR/$SAMPLE_ID.sam.gz
BAM_FILE=$DATA_DIR/$SAMPLE_ID.bam
SORTED_BAM_FILE=$DATA_DIR/$SAMPLE_ID.sort.bam

samtools view -b -o $BAM_FILE $SAM_FILE


if [ "$?" -ne "0" ]; then
    echo "samtools view failed";
    exit 1
fi


samtools sort -T /tmp/$SAMPLE_ID.sorted -o $SORTED_BAM_FILE $BAM_FILE

if [ "$?" -ne "0" ]; then
    echo "samtools sort failed";
    exit 1
fi


samtools index $SORTED_BAM_FILE

if [ "$?" -ne "0" ]; then
    echo "samtools index failed";
    exit 1
fi


rm $BAM_FILE
rm $SAM_FILE

