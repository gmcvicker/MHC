import sys
import gff
import argparse
import chromosome

import util.info

import genome.db


DEFAULT_CHROM_PATH="/home/gm114/data/ucsc/hg19/chromInfo.txt.gz"
DEFAULT_EXON_FILE="/home/gm114/data/MHC/mhc_all_exons.txt"

def parse_options():
    parser = argparse.ArgumentParser()

    parser.add_argument("--chrom", default="6")
    parser.add_argument("--start", default=28000000, type=int,
                        help="start position of provided genotype file")
    
    parser.add_argument("--end", default=34000000, type=int,
                        help="end position of provided genotype file")

    parser.add_argument("--exon_file", default=DEFAULT_EXON_FILE)
    parser.add_argument("--chrom_file", default=DEFAULT_CHROM_PATH)
    parser.add_argument("--min_size", default=0, type=int)

    #    parser.add_argument("--output_dir", required=True)
    
    options = parser.parse_args()

    return options



def main():
    options = parse_options()

    util.info.write_info(sys.stdout, options)

    gdb = genome.db.GenomeDB(assembly="hg19")
    seq_track = gdb.open_track("seq")
    
    chrom_dict = chromosome.get_chromosome_dict(options.chrom_file)

    f = open(options.exon_file)

    for line in f:
        if line.startswith("#"):
            continue
        
        words = line.split()
        
        gene_id = words[0]
        gene_name = words[1]
        chrom_name = words[2]
        exon_num = int(words[3])
        start = int(words[4])
        end = int(words[5])
        strand = int(words[6])

        exon_len = end - start + 1

        if exon_len < options.min_size:
            start = start - (options.min_size - exon_len)/2
            end = end + (options.min_size - exon_len)/2

            sys.stderr.write("extended exon from %d bp to %d bp\n" %
                             (exon_len, end - start + 1))

        
        seq_str = seq_track.get_seq_str(chrom_name, start, end)

        sys.stdout.write(">%s exon %d\n%s\n" % (gene_name, exon_num, seq_str))


main()
            
