#!/bin/bash
#BSUB -n 1
#BSUB -W 0:30
#BSUB -J jobArray[1-60]
#BSUB -o /home/gm114/lsf/%J.%I.out
#BSUB -e /home/gm114/lsf/%J.%I.err
#BSUB -q short


SEQ_DIR=$HOME/data/SGDP/mhc_exon_seqs

SEQ_LIST=$SEQ_DIR/exon_seq_list.txt

SEQ_ID=`sed -n ${LSB_JOBINDEX}p $SEQ_LIST | awk '{print $1}'`

KALIGN=$HOME/kalign/kalign

gunzip -c $SEQ_DIR/$SEQ_ID.fa.gz | $KALIGN | gzip > $SEQ_DIR/$SEQ_ID.aln.gz


if [ "$?" -ne "0" ]; then
    echo "kalign failed";
    exit 1
fi


echo "done"

