#!/bin/sh
#
#
# preserve shell environment vars:
#$ -V
#
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# job needs 8GB of mem:
#$ -l h_vmem=8g

#
# This sript downloads great ape genome project SRA files
# one at a time using the aspera connect client.
# The SRA files must then be converted to fastq files and deleted
# using another script.
#

DIR=/mnt/gluster/data/external_public_supp/great_ape_genome_project/

# ACC_FILE=$DIR/to_download.txt
ACC_FILE=$DIR/SraAccList.txt

while read LINE; do
    ACC=`echo $LINE | awk '{print $1}'`
    
# ACC=`sed -n ${SGE_TASK_ID}p $ACC_FILE | awk '{print $1}'`
    echo "accession: $ACC"


    PFX3=${ACC:0:3}
    PFX6=${ACC:0:6}

    NCBI_PATH=/sra/sra-instant/reads/ByRun/sra/$PFX3/$PFX6/$ACC/$ACC.sra

    ASP_DIR=$HOME/aspera/connect

    $ASP_DIR/bin/ascp -k 1 -T -r -l 300m -i $ASP_DIR/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:$NCBI_PATH $DIR

    if [ "$?" -ne "0" ]; then
	echo "ascp failed for accession $ACC";
	exit 1
    fi
    
    echo "done"
done <$ACC_FILE

