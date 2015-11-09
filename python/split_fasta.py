import sys
import os
import gzip

if len(sys.argv) != 3:
    sys.stderr.write("usage: %s FASTA_FILE OUTPUT_DIR\n")
    exit(2)

fasta_file = sys.argv[1]
output_dir = sys.argv[2]

if fasta_file.endswith(".gz"):
    f = gzip.open(fasta_file, "rb")
else:
    f = open(fasta_file, "r")

os.makedirs(output_dir)

already_seen = set([])

out_f = None

for line in f:
    if line.startswith(">"):
        if out_f:
            out_f.close()

        seq_name = line[1:].rstrip()
        if seq_name in already_seen:
            raise ValueError("duplicate sequence name '%s'\n" % seq_name)
        already_seen.add(seq_name)

        if fasta_file.endswith(".gz"):
            out_filename = "%s/%s.fasta.gz" % (output_dir, seq_name)
            out_f = gzip.open(out_filename, "wb")
        else:
            out_filename = "%s/%s.fasta" % (output_dir, seq_name)
            out_f = open(out_filename, "w")

    if out_f is None:
        raise IOError("first line of file did not start with '>'")
    
    out_f.write(line)

if out_f:
    out_f.close()
        

        

        

        
            


    
