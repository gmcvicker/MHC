import sys
import numpy as np
import tables
import argparse

import util.info

SNP_UNDEF = -1


def parse_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("h5_input_filename")
    
    parser.add_argument("--samples_file",
                        default="/home/gm114/data/SGDP/all_samples.txt")

    parser.add_argument("--region", default=None)
                        
    
    options = parser.parse_args()
    return options


def calc_hetz(gtype_vec):
    defined = gtype_vec != SNP_UNDEF    
    n_defined = 2.0 * np.sum(defined).astype(np.float32)

    if n_defined == 0.0:
        return 0.0
    
    # count number of reference and alternate alleles
    n_alt = np.sum(gtype_vec[defined]).astype(np.float32)
    n_ref = n_defined - n_alt

    # hetz = 2*p*q
    return 2.0 * (n_alt/n_defined) * (n_ref/n_defined)
    


def get_sample_indices(samp_file, region):
    f = open(samp_file)

    header = f.readline()

    index_list = []
    
    for line in f:
        words = line.split("\t")
        samp_region = words[7]

        if samp_region == region:
            samp_idx = int(words[17])
            index_list.append(samp_idx)

    sys.stderr.write("read %d samples from region %s\n" %
                     (len(index_list), region))

    return np.array(index_list)
            

        


        
    
    


def main():
    options = parse_options()

    util.info.write_info(sys.stdout, options)

    if options.region and options.samples_file:
        indices = get_sample_indices(options.samples_file, options.region)
    else:
        indices = None
    
    
    h5f = tables.openFile(options.h5_input_filename, 'r')

    if indices is None:
        gtypes = h5f.root.gtypes[:,:]
    else:
        gtypes = h5f.root.gtypes[:,indices]


    n_sites = gtypes.shape[0]
    n_samp = gtypes.shape[1]

    sys.stderr.write("genotype matrix shape: %s\n" % str(gtypes.shape))

    sys.stderr.write("calculating heterozygosity\n")
    hetz = np.apply_along_axis(calc_hetz, 1, gtypes)
    
    for i in range(n_sites):
        sys.stdout.write("%d %.4f\n" % (i, hetz[i]))
    
    
    h5f.close()


main()
