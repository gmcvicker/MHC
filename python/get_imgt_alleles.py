import sys
import numpy as np
import argparse


class Alignment(object):
    
    def __init__(self):
        self.clear()

    def clear(self):
        self.n_seq = 0
        self.seq_names = []
        self.seqs = []
        self.seq_len = 0

        self.allele2idx = {}
    
        
    def _skip_comments(self, handle):
        # first 5 lines of file are comments
        # skip over first 5 lines, which are comments
        line_num = 0
        for line in handle:
            line_num += 1

            if line_num == 5:
                break
    
            
    def _read_blocks(self, f):
        in_header = False
        i = 0

        for line in f:
            line = line.strip()

            if line == "":
                # start of new alignment block
                in_header = True
                header_line = 0
                i = 0
                # each alignment block header looks like this:
                #
                # cDNA              76
                # AA codon          -4
                #                   |                
            elif in_header:
                header_line += 1
                if header_line == 1:
                    if line.startswith("Please"):
                        # end of file
                        break
                    elif line.startswith("cDNA"):
                        pass
                    else:
                        raise ValueError("expected line to start with"
                                         "'cDNA':\n'%s'\n" % line)
                    # sys.stderr.write("%s\n" % line)
                elif header_line == 2:
                    if not line.startswith("AA"):
                        raise ValueError("expected line to start with 'AA'")
                    # sys.stderr.write("%s\n" % line)
                elif header_line == 3:
                    # last line of header
                    in_header = False
            else:
                # we are in alignment block
                words = line.split()

                # get name of allele, and sequence
                allele_name = words[0]
                sequence = list("".join(words[1:]))

                if i >= self.n_seq:
                    # we are in first alignment block,
                    # add a new allele / sequence
                    self.n_seq += 1

                    self.seq_names.append(allele_name)
                    self.seqs.append(sequence)

                    if allele_name in self.allele2idx:
                        raise ValueError("duplicated allele name %s" % 
                                         allele_name)
                    self.allele2idx[allele_name] = i
                else:
                    # extend existing sequence
                    if self.seq_names[i] != allele_name:
                        # some of the loci have one or two 'longer' alleles
                        # that have extra lines at the end of the file
                        sys.stderr.write("WARNING: expected allele name to "
                                         "be '%s' but got '%s', skipping\n" %
                                         (self.seq_names[i], allele_name))
                    else:
                        self.seqs[i].extend(sequence)
                i += 1

        self.seq_len = len(self.seqs[i])
        
        
    def _set_nucleotides(self):
        for i in range(1, self.n_seq):
            for j in range(self.seq_len):
                if self.seqs[i][j] == "-":
                    # - indicates matches reference
                    self.seqs[i][j] = self.seqs[0][j]
    

                
    def read_alignment(self, filename):
        self.clear()

        f = open(filename)

        self._skip_comments(f)
        self._read_blocks(f)
        sys.stderr.write("%d x %d\n" % (self.n_seq, self.seq_len))
        self._set_nucleotides()


    def write(self, handle, nuc_sep=""):
        for i in range(self.n_seq):
            handle.write("%s %s\n" % (self.seq_names[i], 
                                      nuc_sep.join( self.seqs[i])))

            
    def split_by_exon(self):
        # exons are indicated by '|'
        starts = [0]
        ends = []
        for i in range(self.seq_len):
            if self.seqs[0][i] == '|':
                starts.append(i+1)
                ends.append(i-1)
        ends.append(i)

        exons = []
        for (start, end) in zip(starts, ends):
            ex = Alignment()

            ex.n_seq = self.n_seq
            ex.seq_len = end - start + 1
            ex.seq_names = list(self.seq_names)

            # copy alignment portion of this exon
            for i in range(self.n_seq):
                ex.seqs.append(list(self.seqs[i][start:end+1]))

            exons.append(ex)

        return exons
        
                



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

aln = Alignment()
aln.read_alignment(filename)
# aln.write(sys.stdout)
exons = aln.split_by_exon()

write_exons(options, exons)



