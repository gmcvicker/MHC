
import sys
import gff
import argparse
import chromosome

import util.info

DEFAULT_GFF_PATH="/home/gm114/data/GENCODE/hg19/gencode.v19.annotation.gff3.gz"
DEFAULT_CHROM_PATH="/home/gm114/data/ucsc/hg19/chromInfo.txt.gz"

def parse_options():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--chrom", default="chr6")
    parser.add_argument("--start", default=28000000, type=int)
    parser.add_argument("--end", default=34000000, type=int)
    parser.add_argument("--gff", default=DEFAULT_GFF_PATH)
    parser.add_argument("--chrom_file", default=DEFAULT_CHROM_PATH)

    # parser.add_argument("output_filename")

    options = parser.parse_args()

    return options



def main():
    options = parse_options()

    util.info.write_info(sys.stdout, options)
    
    chrom_dict = chromosome.get_chromosome_dict(options.chrom_file)
    
    gene_dict, tr_dict, gene_chrom_dict, tr_chrom_dict = \
      gff.read_gff(options.gff, chrom_dict,
                   region_chrom=options.chrom,
                   region_start=options.start,
                   region_end=options.end)

    for gene in gene_chrom_dict[options.chrom]:
        exons = gene.get_merged_exons()

        i = 0

        if gene.strand == -1:
            exons  = exons[::-1]
        
        for ex in exons:
            i += 1
            sys.stdout.write("%s %s %s %d %d %d %d\n" % 
                             (gene.gene_id, gene.gene_name, gene.chrom.name, i,
                              ex.start, ex.end, gene.strand))

main()
            
