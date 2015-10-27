import sys
import gzip


MHC_CONTIGS=set(["chr6_GL000250v2_alt", "chr6_GL000251v2_alt",
                 "chr6_GL000252v2_alt", "chr6_GL000253v2_alt",
                 "chr6_GL000254v2_alt", "chr6_GL000255v2_alt",
                 "chr6_GL000256v2_alt"])

MHC_START = 29000000
MHC_END = 34000000

class Alignment(object):
    def __init__(self, line):
        words = line.split("\t")
        self.query_id = words[0]
        self.subject_id = words[1]
        self.perc_identity = float(words[2])
        self.aln_len = int(words[3])
        self.mismatch_count = int(words[4])
        self.gap_open_count = int(words[5])
        self.query_start = int(words[6])
        self.query_end = int(words[7])
        self.subject_start = int(words[8])
        self.subject_end = int(words[9])
        self.e_val = float(words[10])
        self.bit_score = float(words[11])


    def is_mhc(self):
        # return true if alignment subject is MHC region
        if(self.subject_id == "chr6" and
           self.subject_end > MHC_START and
           self.subject_start < MHC_END):
            # is aligned to MHC region of chr6
            return True
        elif self.subject_id in MHC_CONTIGS:
            # is aligned to one of the MHC alt contigs in assembly...
            return True

        # aligned to other region...
        return False

    def subject_str(self):
        return "%s:%d-%d" % (self.subject_id, self.subject_start, self.subject_end)
        
                

class Contig(object):
    def __init__(self, locus, short_samp_id, name):
        self.locus = locus
        self.short_samp_id = short_samp_id
        self.name = name      
        self.length = int(name.split("_")[3])
        self.best_mhc_align = None
        self.best_non_mhc_align = None
        self.n_mhc_align = 0
        self.n_non_mhc_align = 0


    def write(self, f):
        """writes information about this contig to provided file handle"""
        if self.best_mhc_align:
            mhc_align_str = self.best_mhc_align.subject_str()
            mhc_score_str = str(self.best_mhc_align.bit_score)
        else:
            mhc_align_str = "."
            mhc_score_str = "0"

        if self.best_non_mhc_align:
            non_mhc_align_str = self.best_non_mhc_align.subject_str()
            non_mhc_score_str = str(self.best_non_mhc_align.bit_score)
        else:
            non_mhc_align_str = "."
            non_mhc_score_str = "0"
                                                          
        f.write("\t".join([self.locus, self.short_samp_id, self.name,
                           str(self.length), mhc_align_str, non_mhc_align_str,
                           mhc_score_str, non_mhc_score_str,
                           str(self.n_mhc_align), str(self.n_non_mhc_align)]) + "\n")
    


def main():

    if len(sys.argv) < 2:
        sys.stderr.write("usage: %s contig_blast_output1 "
                         "[contig_blast_output2 ...]\n" % sys.argv[0])
        exit(2)

    sys.stdout.write("\t".join(["LOCUS", "SHORT.ID", "CONTIG.NAME", "CONTIG.LENGTH",
                                "MHC.BEST.HIT", "NON.MHC.BEST.HIT",
                                "MHC.BEST.SCORE", "NON.MHC.BEST.SCORE",
                                "N.MHC.HIT", "N.NON.MHC.HIT"]) + "\n")
        
    for path in sys.argv[1:]:
        words = path.split("/")
        filename = words[-1]
        parent_dir = words[-2]

        # parent directory name looks like:
        # "928-DPB1_S25_L001"
        # this is SHORT.SAMPLE.ID followed by LOCUS followed by
        # illumina samp/lane number
        identifier =  parent_dir.split("_", )[0]
        if "-" in identifier:
            short_samp_id, locus = identifier.split("-", 1)
        else:
            short_samp_id = "."
            locus = identifier

        f = gzip.open(path)

        cur_contig = None

        for line in f:
            align = Alignment(line)

            if cur_contig is None:
                cur_contig = Contig(locus, short_samp_id, align.query_id)
            elif cur_contig.name != align.query_id:
                # start of alignments for new contig, write out current one
                cur_contig.write(sys.stdout)
                cur_contig = Contig(locus, short_samp_id, align.query_id)

            if align.is_mhc():
                # this is an MHC alignment
                cur_contig.n_mhc_align += 1
                if (cur_contig.best_mhc_align is None or
                    cur_contig.best_mhc_align.bit_score < align.bit_score):
                    # scores better than previously seen alignments to MHC
                    cur_contig.best_mhc_align = align
            else:
                # this is alignment to non-MHC part of genome
                cur_contig.n_non_mhc_align += 1
                if (cur_contig.best_non_mhc_align is None or
                    cur_contig.best_non_mhc_align.bit_score < align.bit_score):
                    cur_contig.best_non_mhc_align = align

        if cur_contig:
            cur_contig.write(sys.stdout)
        
        f.close()
        

        
main()        
        
        
        
        
