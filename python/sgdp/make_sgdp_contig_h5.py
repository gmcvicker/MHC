import sys
import tables

import pysam



SAMPLES_FILE = "/home/gm114/data/SGDP/C-team-mag/261.samples"

CONTIG_H5 = "/home/gm114/data/SGDP/mag_mhc_contigs.h5"
DEPTH_H5 = "/home/gm114/data/SGDP/mag_mhc_depth.h5"

MHC_START = 29000000
MHC_END = 34000000

N_SAMPLES = 261



def create_contig_h5_file():
    # initialize HDF5 output file
    h5f = tables.openFile(CONTIG_H5, "w")
    atom = tables.UInt8Atom(dflt=ord('.'))

    n_sites = MHC_END - MHC_START + 1
    shape = (n_sites, N_SAMPLES)
    
    carray = h5f.createCArray(h5f.root, "contigs", atom, shape,
                              filters=tables.Filters(complevel=1, complib="zlib"))


    return carray



def create_depth_h5_file():
    h5f = tables.openFile(DEPTH_H5, "w")
    atom = tables.Int8Atom(dflt=0)

    n_sites = MHC_END - MHC_START + 1
    shape = (n_sites, N_SAMPLES)

    carray = h5f.createCArray(h5f.root, "depth", atom, shape,
                              filters=tables.Filters(complevel=1, complib="zlib"))


    
    


def main():
    f = open(SAMPLES_FILE, "r")

    contig_h5 = create_contig_h5_file()
    depth_h5 = create_depth_h5_file()

    depth_carray = depth_h5.getNode("/depth")


    mhc_len = MHC_END - MHC_START + 1

    indv_index = 0
    
    for l in f:
        words = l.rstrip().split()
        cteam_id = words[1]
        
        bam_filename = "/home/gm114/data/SGDP/mag_mhc_mapped/%s.sort.bam" % cteam_id
        bamfile = pysam.Samfile(bam_filename, "rb")

        sys.stderr.write("%s\n" % cteam_id)
        depth_array = np.zeros(mhc_len)
        
        for read in bamfile.fetch("chr6", MHC_START, MHC_END):
            for blk in read.get_blocks():
                start = blk[0] - MHC_START
                end = blk[1] - MHC_START

                if start < 0:
                    start = 0
                if end < 0:
                    end = 0

                if start > mhc_len:
                    start = mhc_len
                if end > mhc_len:
                    end = mhc_len

                depth_array[start:end] = 
                    
                print "%d-%d" % (blk[0], blk[1])

        indv_index = 0

                
                
            
            
        

        bamfile.close()

    
    contig_h5.close()
    depth_h5.close()
        
    f.close()
    


main()
