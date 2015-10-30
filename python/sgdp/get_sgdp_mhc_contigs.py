import sys
import pysam


MHC_CHROM = "chr6"
MHC_START = 2800000
MHC_END = 34000000

PATCH_FILE = "/home/gm114/data/GRCh38/Patches_and_alternate_loci_for_GRCh38.p1.tsv"

def read_alt_loci():
    f = open(PATCH_FILE)

    mhc_alt_loci = set([])
    
    for line in f:
        words = line.rstrip().split()

        alt_type = words[2]
        if alt_type == "MHC":
            name = "chr6_" + words[0].split("|")[0].replace(".", "v") + "_alt"
            mhc_alt_loci.add(name)

    return mhc_alt_loci


            



def main():
    if len(sys.argv) != 3:
        sys.stderr.write("usage: %s <indexed bamfile> <sample_id> | gzip > "
                         "<contigs.fastq.gz>\n" % sys.argv[0])
        exit(2)

    mhc_alt_loci = read_alt_loci()
    
    filename = sys.argv[1]
    sample_id = sys.argv[2]
    
    bamfile = pysam.Samfile(filename, "rb")

    mhc_ref_names = set([])
    
    for ref in bamfile.references:
        if ref in mhc_alt_loci or ref.startswith("HLA"):
            mhc_ref_names.add(ref)

    for ref in mhc_ref_names:
        for read in bamfile.fetch(ref):
            # pysam repeats reads, possibly for printing PE reads,
            # should fix this so only written once
            sys.stdout.write("@%s %s\n%s\n+\n%s\n" % (read.qname, sample_id, 
                                                      read.query, 
                                                      read.qqual))
                

            
    bamfile.close()
    


main()
