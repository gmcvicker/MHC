
FILE=/home/gm114/data/SGDP/VCFS/LP6005442-DNA_G11.6.filtered.vcf.gz

MHC_START=28000000
MHC_END=34000000

gunzip -c $FILE | egrep -v "^#" | awk -v start=$MHC_START -v end=$MHC_END 'BEGIN{print start, end} ($2 > start && $2 < end) {print $1, $2, $10}'



