#!/bin/sh
#
#
# preserve shell environment vars:
#$ -V
#
# run a job-array consisting of tasks 1-399
#$ -t 1-399
#
# write output and err files to current dir
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# jobs need 2GB of mem and lots of IO:
#$ -l h_vmem=2g,bigio=2


# Ths script converts sra files to fastq files
# WARNING: This script REMOVES the sra file to save space

DIR=/mnt/gluster/data/external_public_supp/great_ape_genome_project/

# ACC_FILE=$DIR/to_download.txt
ACC_FILE=$DIR/SraAccList.txt

ACC=`sed -n ${SGE_TASK_ID}p $ACC_FILE | awk '{print $1}'`

SRA_DIR=$HOME/sratoolkit.2.5.2-ubuntu64

echo $ACC

$SRA_DIR/bin/fastq-dump -I --outdir $DIR --split-files --gzip $DIR/$ACC.sra

if [ "$?" -ne "0" ]; then
    echo "fastq-dump failed";
    exit 1
fi

# REMOVE FILE AFTER EXTRACTING FASTQ
rm $DIR/$ACC.sra

echo "done"

