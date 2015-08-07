#!/bin/bash
#BSUB -n 1
#BSUB -W 0:10
#BSUB -J jobArray[1-261]
#BSUB -o /home/gm114/lsf/%J.%I.out
#BSUB -e /home/gm114/lsf/%J.%I.err
#BSUB -q short

CONTIG_DIR=/n/data1/hms/genetics/reich/1000Genomes/lh3/C-team-mag

SAMPLES_FILE=$CONTIG_DIR/261.samples

SAMPLE_ID=`sed -n ${LSB_JOBINDEX}p $SAMPLES_FILE | awk '{print $2}'`

echo $SAMPLE_ID >&2

BAM_FILE=$CONTIG_DIR/hs38/$SAMPLE_ID.srt.bam

OUT_DIR=/home/gm114/data/SGDP/mag_mhc_mapped
OUT_FILE=$OUT_DIR/$SAMPLE_ID.fq.gz

SCRIPT=/home/gm114/proj/MHC/python/get_sgdp_mhc_contigs.py

python $SCRIPT $BAM_FILE $SAMPLE_ID | gzip > $OUT_FILE

if [ "$?" -ne "0" ]; then
    echo "$SCRIPT failed" >&2
    exit 1
fi


