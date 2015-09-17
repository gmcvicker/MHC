import sys
import argparse
import pysam
import os
import util.info
import gzip


DEFAULT_EXON_FILE = "/home/gm114/data/GENCODE/hg38/mhc_all_exons.txt"
DEFAULT_CLASSIC_GENE_FILE = "/home/gm114/data/GENCODE/hg38/mhc_classical_genes.txt"
BAM_DIR = "/home/gm114/data/SGDP/mag_mhc_mapped"



class Exon(object):

    def __init__(self, chrom_name, start, end, strand=0):
        self.chrom = chrom_name
        self.start = start
        self.end = end
        self.strand = strand

        self.gene_name = None
        self.gene_id = None
        self.num = 0

    
        

def parse_options():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--exon_file", default=DEFAULT_EXON_FILE)
    parser.add_argument("--gene_subset", default=DEFAULT_CLASSIC_GENE_FILE)
    parser.add_argument("--bam_dir", default=BAM_DIR)

    parser.add_argument("output_dir")

    options = parser.parse_args()

    return options



def read_exons(options):
    
    f = open(options.gene_subset)
    gene_subset = set([])
    for line in f:
        gene_subset.add(line.strip())
    f.close()    

    f = open(options.exon_file)

    exon_list = []
    
    for line in f:
        if line.startswith("#"):
            continue

        words = line.strip().split()
        gene_id = words[0]
        gene_name = words[1]
        chrom_name = words[2]
        exon_num = int(words[3])
        start = int(words[4])
        end = int(words[5])
        strand = int(words[6])

        if gene_name in gene_subset:
            exon = Exon(chrom_name, start, end, strand=strand)
            exon.gene_name = gene_name
            exon.gene_id = gene_id
            exon.num = exon_num

            exon_list.append(exon)

    f.close()

    return exon_list

    


def main():
    options = parse_options()
    
    info_f = open(options.output_dir + "/info.txt", "w")
    util.info.write_info(info_f, options)
    info_f.close()

    exons = read_exons(options)

    bam_filenames = []
    for filename in os.listdir(options.bam_dir):
        if filename.endswith(".sort.bam"):
            bam_filenames.append(filename)

    for exon in exons:
        exon_str = "%s.%d" % (exon.gene_name, exon.num)
        sys.stderr.write("%s\n" % exon_str)
        
        output_filename = options.output_dir + "/%s.fa.gz" % exon_str
        out_f = gzip.open(output_filename, "wb")
        
        for bam_filename in bam_filenames:
            bam_f = pysam.Samfile(options.bam_dir + "/" + bam_filename, "rb")

            sample_id = bam_filename.split(".")[0]

            last = ""
            for read in bam_f.fetch(exon.chrom, exon.start, exon.end):
                if read.qname != last:
                    out_f.write(">%s %s\n%s\n" % (read.qname, sample_id, 
                                                  read.query))
                    last = read.qname

        out_f.close()

                    

main()
    
