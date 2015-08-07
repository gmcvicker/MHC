
import sys
import gff
import argparse
import chromosome

import util.info


# DEFAULT_GFF_PATH="/home/gm114/data/GENCODE/hg19/gencode.v19.annotation.gff3.gz"
# DEFAULT_CHROM_PATH="/home/gm114/data/ucsc/hg19/chromInfo.txt.gz"

DEFAULT_GFF_PATH="/home/gm114/data/GENCODE/hg38/gencode.v23.annotation.gff3.gz"
DEFAULT_CHROM_PATH="/home/gm114/data/ucsc/hg38/chromInfo.txt.gz"


def parse_options():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--chrom", default="chr6")
    parser.add_argument("--start", default=28000000, type=int)
    parser.add_argument("--end", default=34000000, type=int)
    parser.add_argument("--gff", default=DEFAULT_GFF_PATH)
    parser.add_argument("--chrom_file", default=DEFAULT_CHROM_PATH)

    parser.add_argument("exon_output_filename")
    parser.add_argument("gene_output_filename")

    options = parser.parse_args()

    return options



def main():
    options = parse_options()

    exon_out_f = open(options.exon_output_filename, "w")
    gene_out_f = open(options.gene_output_filename, "w")
    
    util.info.write_info(exon_out_f, options)
    util.info.write_info(gene_out_f, options)
        
    chrom_dict = chromosome.get_chromosome_dict(options.chrom_file)
    
    gene_dict, tr_dict, gene_chrom_dict, tr_chrom_dict = \
      gff.read_gff(options.gff, chrom_dict,
                   region_chrom=options.chrom,
                   region_start=options.start,
                   region_end=options.end)

    gene_num = 0
    for gene in gene_chrom_dict[options.chrom]:
        exons = gene.get_merged_exons()

        gene_num += 1
        gene_out_f.write("%s %s %s %d %d %d %d\n" % (gene.gene_id, gene.gene_name, 
                                                     gene.chrom.name, gene_num, gene.start,
                                                     gene.end, gene.strand))
        exon_num = 0

        if gene.strand == -1:
            exons = exons[::-1]
        
        for ex in exons:
            exon_num += 1
            exon_out_f.write("%s %s %s %d %d %d %d\n" % 
                             (gene.gene_id, gene.gene_name, gene.chrom.name, exon_num,
                              ex.start, ex.end, gene.strand))


    exon_out_f.close()
    gene_out_f.close()

main()
            
