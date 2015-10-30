import sys
import gzip

import os
import re

import argparse

import pysam
import numpy as np

import tables

SNP_UNDEF = -1

HOME = os.environ['HOME']
DATA_DIR = HOME + "/data/great_ape_genome_project"





def read_samples(species):
    sample_file = DATA_DIR + "/samples.%s.txt" % species

    sample_list = []
    f = open(sample_file)
    for line in f:
        sample_list.append(line.rstrip())
            
    return sample_list




def get_samp_genos(filename, n_samp, chrom, start, end):
    n_sites = end - start + 1

    # unfortunately VCFs provided by Ape genome project do not
    # indicate which sites are genotyped... assume that
    # if no genotype provided that matches reference...
    genos = np.empty((n_sites, n_samp), dtype=np.int8)
    genos[:] = 0

    f = gzip.open(filename)

    for line in f:
        if line.startswith("#"):
            continue

        words = line.split()

        chrom_str = words[0]
        pos = int(words[1])
        ref_base = words[2]
        alt_base = words[3]

        if chrom_str != chrom:
            raise ValueError("script assumes only chr6 is in VCF '%s' != '%s'" %
                             (chrom, chrom_str))

        if pos < start or pos > end:
            raise ValueError("script assumes only MHC region is in VCF")

        vals = words[9:]

        if len(vals) != n_samp:
            raise ValueError("expected vals for %d samples, got %d" %
                             (n_samp, len(vals)))
        
        pos_idx = pos - start

        for samp_idx in range(n_samp):
            gtype_str = vals[samp_idx][0:3]
            if gtype_str[1] != '/':
                raise ValueError("expected genotype string like 0/0, 0/1, or 1/1")
            
            # count number of non-reference alleles...
            if '.' in gtype_str:
                genos[pos_idx, samp_idx] = SNP_UNDEF
            else:
                gt1 = int(gtype_str[0])
                gt2 = int(gtype_str[2])

            nonref_count = 0

            if gt1 > 0:
                nonref_count += 1
            if gt2 > 0:
                nonref_count += 1
            
            genos[pos_idx, samp_idx] = nonref_count


    f.close()
    sys.stderr.write("\n")

    return genos





def parse_options():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--chrom", default="chr6")
    parser.add_argument("--start", default=28000000)
    parser.add_argument("--end", default=34000000)

    parser.add_argument("output_dir")

    options = parser.parse_args()

    return options
    
        
    
def main():
    options = parse_options()
    
    species_path = DATA_DIR + "/species.txt"

    species_list = []
    f = open(species_path)
    for l in f:
        species_list.append(l.rstrip())
    
    for species in species_list:
        sys.stderr.write("%s\n" % species)    
        
        samples = read_samples(species)
        
        n_samp = len(samples)
        n_sites = options.end - options.start + 1

        output_filename = options.output_dir + "/%s_mhc_genotypes.h5" % species
        sys.stderr.write("writing to file %s\n" % output_filename)

        # gtypes = np.zeros((n_sites, n_samp), dtype=np.int8)

        # initialize HDF5 output file
        h5f = tables.openFile(output_filename, "w")
        atom = tables.Int8Atom(dflt=SNP_UNDEF)
        shape = (n_sites, n_samp)
        carray = h5f.createCArray(h5f.root, "gtypes", atom, shape,
                                  filters=tables.Filters(complevel=1, complib="zlib"))

        vcf_path = DATA_DIR + "/VCFs/SNPs/%s.chr6_mhc.vcf.gz" % species
        
        carray[:] = get_samp_genos(vcf_path, n_samp, options.chrom, 
                                   options.start, options.end)
        
        h5f.flush()
        h5f.close()

      
main()
