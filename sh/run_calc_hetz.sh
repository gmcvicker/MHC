#!/bin/bash
#BSUB -n 1
#BSUB -W 8:00
#BSUB -J jobArray[1-2]
#BSUB -o /home/gm114/lsf/%J.%I.out
#BSUB -e /home/gm114/lsf/%J.%I.err
#BSUB -q short


for REGION in Africa America CentralAsiaSiberia EastAsia Oceania SouthAsia WestEurasia; do
    echo $REGION
    python $HOME/proj/MHC/python/calc_hetz.py --region $REGION /home/gm114/data/SGDP/mhc_genotypes.h5 | gzip > mhc_hetz.$REGION.txt.gz
done


	      

