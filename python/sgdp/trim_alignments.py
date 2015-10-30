import sys
import gzip
import os
import numpy as np
import util.info


ALN_DIR="/home/gm114/data/SGDP/mhc_exon_seqs"


def read_seqs(path):
    cur_seq = ""
    seq_list = []
    header_list = []

    f = gzip.open(path)

    for line in f:
        line = line.strip()
        if line.startswith(">"):
            header_list.append(line)
            if cur_seq:
                seq_list.append([ord(x) for x in cur_seq])
            cur_seq = ""
        else:
            cur_seq = cur_seq + line

    if cur_seq:
        seq_list.append([ord(x) for x in cur_seq])


    return header_list, np.array(seq_list, dtype=np.uint8)



def count_gaps(seq_array):
    # count number of gaps at each column of alignment

    gap_val = ord("-")

    return np.apply_along_axis(np.sum, 0, seq_array == gap_val)

    
    # n_col = len(seq_list[0])
    # n_row = len(seq_list)
    
    # gap_count = np.zeros(n_col, dtype=np.int16)

    # sys.stderr.write("  %d rows %d cols\n" % (n_row, n_col)
    # for cur_seq in seq_list:
    #     for i in range(len(cur_seq)):
    #         if cur_seq[i] == "-":
    #             gap_count[i] += 1
    # return gap_count
                
                

        
def main():        

    util.info.write_info(sys.stderr)
    
    for filename in os.listdir(ALN_DIR):
        if filename.endswith(".aln.gz"):
            sys.stderr.write(filename + "\n")
            path = ALN_DIR + "/" + filename

            sys.stderr.write("  reading alignment\n")
            header_list, seq_array = read_seqs(path)

            sys.stderr.write("seq_array.shape: " + str(seq_array.shape))

            sys.stderr.write("  counting gaps\n")
            gap_counts = count_gaps(seq_array)
            # sys.stderr.write(" ".join([str(x) for x in gap_counts]) + "\n")

            # find portions of the alignment that have low number of gaps
            n_row = seq_array.shape[0]
            n_col = seq_array.shape[1]
            gap_frac = gap_counts.astype(np.float32) / n_row

            smooth_win_sz = 5
            smooth_win = np.array([1.0] * smooth_win_sz)
            smooth_gap_frac = np.convolve(gap_frac, smooth_win, 'same') / smooth_win_sz
            
            low_gap_idx = np.nonzero(smooth_gap_frac < 0.25)[0]

            # trim alignment to span first and last low-gap portions of alignment
            if low_gap_idx.shape[0] > 0:
                low_gap_start = low_gap_idx[0]
                low_gap_end = low_gap_idx[-1]

                aln_start = max(0, low_gap_start - smooth_win_sz)
                aln_end   = min(n_col, low_gap_end + smooth_win_sz)

                for i in range(n_row):
                    sys.stderr.write("".join([chr(x) for x in seq_array[i, aln_start:aln_end+1]]) + "\n")
                
            
            #sys.stderr.write(" ".join(["%.3f" % x for x in smooth_gap_frac]) + "\n")
            #sys.stderr.write(" ".join(["%d" % x for x in low_gap]) + "\n")

            

            
    
main()
