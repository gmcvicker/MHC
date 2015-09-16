#!/bin/bash
#BSUB -n 1
#BSUB -W 12:00
#BSUB -J jobArray[1-396]%5
#BSUB -o /home/gm114/lsf/%J.%I.out
#BSUB -e /home/gm114/lsf/%J.%I.err
#BSUB -q short

DIR=$HOME/data/great_ape_genome_project/fastq

ACC_FILE=$DIR/to_download.txt

acc=`sed -n ${LSB_JOBINDEX}p $ACC_FILE | awk '{print $1}'`

echo $acc
fastq-dump -I --outdir $DIR --split-files --gzip $acc


if [ "$?" -ne "0" ]; then
    echo "fastq-dump failed";
    exit 1
fi


echo "done"


