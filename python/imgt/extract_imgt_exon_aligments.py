import sys
import numpy as np
import argparse

import imgt
                



def parse_options():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("output_dir")
    parser.add_argument("--gene", default="DRB")

    options = parser.parse_args()

    return options
    


def write_exons(options, exons):

    exon_number = 0
    
    for ex in exons:
        exon_number += 1

        out_filename = "%s/%s_exon%d.txt" % \
          (options.output_dir, options.gene, exon_number)

        sys.stderr.write("writing to file %s\n" % out_filename)

        out_f = open(out_filename, "w")
        ex.write(out_f, nuc_sep=" ")
        out_f.close()
    
    

    

def make_diff_matrix(aln):
    seq_array = np.array(aln.seqs)
    
    ##### Currently not used
    diff_matrix = np.zeros((aln.n_seq, aln.n_seq))
    for i in range(aln.n_seq):
        for j in range(i):
            # don't want untyped sites, which are flagged with '*'
            subset = (seq_array[i,] != '*') | (seq_array[j,] != '*')

            diff_matrix[i,j] = np.sum(seq_array[i,subset] != seq_array[j,])
            diff_matrix[j,i] = diff_matrix[i,j]



options = parse_options()

filename = "/home/gm114/data/IMGT/alignments/" + options.gene + "_nuc.txt"

aln = imgt.Alignment()
aln.read_alignment(filename)
# aln.write(sys.stdout)
exons = aln.split_by_exon()

write_exons(options, exons)



