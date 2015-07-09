import sys
import numpy as np
import tables
import argparse

SNP_UNDEF = -1


def parse_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("h5_input_filename")
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
    


def main():
    options = parse_options()
    
    h5f = tables.openFile(options.h5_input_filename, 'r')
    
    gtypes = h5f.root.gtypes[:,:]

    n_sites = gtypes.shape[0]
    n_samp = gtypes.shape[1]

    sys.stderr.write("genotype matrix shape: %s\n" % str(gtypes.shape))

    sys.stderr.write("calculating heterozygosity\n")
    hetz = np.apply_along_axis(calc_hetz, 1, gtypes)
    
    for i in range(n_sites):
        sys.stdout.write("%d %.4f\n" % (i, hetz[i]))
    
    
    h5f.close()


main()
