
import sys
import os
import shutil


assembly_dir = "/home/gm114/data/GRCh38/"
hap_file = assembly_dir + "Patches_and_alternate_loci_for_GRCh38.p1.tsv"
chrom_dir = assembly_dir + "chroms/"

MHC_dir = "/home/gm114/data/GRCh38/MHC/"
os.mkdir(MHC_dir)


f = open(hap_file)

mhc_seq_ids = set([])

# pull out identifiers of MHC sequences
header = f.readline()
for line in f:
    words = line.split()
    seq_id_vers = words[0].split("|")[0]
    gb_acc = words[0].split("|")[1]
    seq_type = words[1]
    region_name = words[2]

    seq_id = seq_id_vers.split(".")[0]
    seq_vers = seq_id_vers.split(".")[1]
    
    if region_name == "MHC":
        sys.stderr.write("  %s %s\n" % (seq_id, gb_acc))
        mhc_seq_ids.add(seq_id)


# find the MHC sequence files and copy them to new dir
for filename in os.listdir(chrom_dir):
    if filename.startswith("chr") and ("_" in filename):
        
        seq_id_vers = filename.split("_")[1]

        seq_id = seq_id_vers.split("v")[0]
        vers = seq_id_vers.split("v")[1]

        if seq_id in mhc_seq_ids:
            sys.stderr.write("%s => %s\n" % (seq_id, filename))
            frm_path = chrom_dir + filename
            to_path = MHC_dir + filename
            shutil.copy(frm_path, to_path)

        

            

        
        
        
        

        
        
    

    
