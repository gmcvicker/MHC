import sys
import gzip

import os
import re

import argparse

import pysam
import numpy as np

import pandas as pd

import tables

SNP_UNDEF = -1

HOME = os.environ['HOME']
DATA_DIR = HOME + "/data/SGDP/VCFS"



class Sample(object):

    def __init__(self, samp_id, name, idx, prefix):

        # identifier string for sample e.g. "LP6005442-DNA_H09"
        self.samp_id = samp_id

        # name of sample, e.g. "S_Bengali-2"
        self.name = name

        # index of sample, first sample is 0
        self.idx = idx

        # prefix of path to datafile for sample, 
        # e.g. /home/gm114/data/SGDP/VCFS/LP6005442-DNA_G11
        self.prefix = prefix




def read_samples():
    sample_file = HOME + "/data/SGDP/samples/hetfa.dblist"
    
    f = open(sample_file)

    idx = 0

    samp_list = []
    
    for line in f:
        cols = line.split()
        samp_name = cols[0]
        samp_id = cols[1]
        samp_path = cols[2]

        words = samp_path.split("/")
        filename = words[-1].replace(".raw.hetfa", "")
        vcf_dir = HOME + "/data/SGDP/samples/" + "/".join(words[-3:-2]) + "/VCFS"

        vcf_prefix = vcf_dir + "/" + filename

        sys.stderr.write(vcf_prefix + "\n")

        samp = Sample(samp_id, samp_name, idx, vcf_prefix)

        samp_list.append(samp)
        idx += 1


    
    return samp_list




def write_output(handle, sample_list, gtype_array):
    samp_ids = " ".join([samp.id for samp in sample_list])

    # write a header line
    handle.write(samp_ids + "\n")

    for row in range(gtype_array.shape[0]):
        handle.write(" ".join(["%d" % x for x in gtype_array[i,:]]) + "\n")

        
        
        
    


def get_samp_genos(filename, chrom, start, end):
    n_sites = end - start + 1
    genos = np.empty(n_sites, dtype=np.int8)
    genos[:] = SNP_UNDEF

    
    try:
        tbx = pysam.Tabixfile(filename)
    except IOError as e:
        sys.stderr.write("WARNING: could not open file %s: %s" % 
                         (filename, e.strerror))
        return genos

    # tabix fetch uses half-open, zero-based indexing
    for row in tbx.fetch(chrom, start-1, end):
        cols = row.split()
        pos = int(cols[1])

        if pos < start or pos > end:
            raise ValueError("position %d is outside of expected range %d-%d" %
                             (pos, start, end))

        # genotype string looks like 0/0:44
        # where 0s indicate ref, 1 (or rarely 2) indicates non-ref
        # for unknown genotype there can be '.'
        # the last number is read depth
        gt1_str = cols[9][0]
        gt2_str = cols[9][2]

        if gt1_str == "." or gt2_str == ".":
            # undefined genotype, leave as -1
            pass
        else:
            gt1 = int(gt1_str)
            gt2 = int(gt2_str)

            nonref_count = 0

            if gt1 > 0:
                nonref_count += 1
            if gt2 > 0:            
                nonref_count += 1

            pos_idx = pos - start
            genos[pos_idx] = nonref_count

            if pos_idx % 100000 == 0:
                sys.stderr.write(".")

    sys.stderr.write("\n")

    tbx.close()

    return genos





def parse_options():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--chrom", default="6")
    parser.add_argument("--start", default=28000000)
    parser.add_argument("--end", default=34000000)
    parser.add_argument("--subset", type=int, default=0)

    parser.add_argument("output_filename")

    options = parser.parse_args()

    return options
    
        
    
def main():
    options = parse_options()
    
    samples = read_samples()

    if (options.subset > 0) and (options.subset < len(samples)):
        # take a subset of the samples
        sys.stderr.write("using subset of %d/%d samples\n" % 
                         (options.subset, len(samples)))
        samples = samples[0:options.subset]
        
    n_samp = len(samples)
    n_sites = options.end - options.start + 1

    # gtypes = np.zeros((n_sites, n_samp), dtype=np.int8)

    # initialize HDF5 output file
    h5f = tables.openFile(options.output_filename, "w")
    atom = tables.Int8Atom(dflt=SNP_UNDEF)
    shape = (n_sites, n_samp)
    carray = h5f.createCArray(h5f.root, "gtypes", atom, shape,
                              filters=tables.Filters(complevel=1, complib="zlib"))
    
    for i in range(n_samp):
        samp = samples[i]

        path = samp.prefix + "." + options.chrom + ".filtered.vcf.gz"
        sys.stderr.write("reading from file %s\n" % path)

        carray[:, i] = \
          get_samp_genos(path, options.chrom, options.start, options.end)
    

    h5f.flush()
    h5f.close()

      
main()
