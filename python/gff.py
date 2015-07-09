
import gzip
import re
import sys

import genome.transcript
import genome.gene
import genome.coord


def read_gff(filename, chrom_dict, strip_version=True,
             region_chrom=None, region_start=None, region_end=None):
    """Reads transcripts and genes from GFF file, returns them as
    dictionaries keyed on gene identifier, transcript identifier, 
    chromosome name. Dictionaries keyed on chromosome name
    contain lists sorted by start coordinate. Chrom, start, 
    end arguments can be used to restrict to specific region.
    """

    # possible feature types in ensembl GFF file are:
    #  CDS
    #  exon
    #  gene
    #  Selenocysteine
    #  start_codon
    #  stop_codon
    #  transcript
    #  UTR

    if filename.endswith(".gz"):
        f = gzip.open(filename)
    else:
        f = open(filename)

    gene_dict = {}
    tr_dict = {}
    tr_ex_dict = {}
    gene_chrom_dict = dict([(x, []) for x in chrom_dict.keys()])
    tr_chrom_dict = dict([(x, []) for x in chrom_dict.keys()])
    

    for line in f:
        if line.startswith("#"):
            # comment line
            continue

        words = line.rstrip().split(None, 8)

        chrom_name = words[0]

        if region_chrom and (chrom_name != region_chrom):
            # skip this gene, not on chrom of interest
            continue
                
        annot_source = words[1]
        feature_type = words[2]

        
        start = int(words[3])
        if region_start and start < region_start:
            # skip feature that is not entirely in specified region
            continue

        end = int(words[4])
        if region_end and end > region_end:
            # skip feature that is not entirely in specified region
            continue
        
        score = words[5]

        strand = 0
        if words[6] == '+':
            strand = 1
        elif words[6] == '-':
            strand = -1

        frame = words[7]

        # last fields are colon-delimited 'attributes' 
        attrib =  words[8].split(";")[:-1]

        attrib_dict = {}
        for a in attrib:
            # may be whitespace or '=' separated
            tok = re.split("=| |\t", a)
            if len(tok) != 2:
                raise ValueError("expected each attribute "
                                 "to have 2 tokens, got: '%s'\n"
                                 "LINE: '%s'" % (a, line))
            
            attrib_dict[tok[0]] = tok[1].replace('"', '')
        
        # sys.stderr.write("%s\n" % feature_type)
        # for k,v in attrib_dict.items():
        #    sys.stderr.write("  %s => %s\n" % (k,v))

        if feature_type == 'transcript':
            # create new transcript
            tr_id = attrib_dict['transcript_id']
            if strip_version:
                tr_id = tr_id.split(".")[0]
            
            tr = genome.transcript.Transcript(name=tr_id)
            tr.chrom = chrom_dict[chrom_name]
            tr.start = start
            tr.end = end
            tr.strand = strand

            if tr_id in tr_dict:
                raise ValueError("Transcript %s defined "
                                 "multiple times in file" % tr_id)

            tr.transcript_id = tr_id

            if strip_version:
                tr.gene_id = attrib_dict['gene_id'].split(".")[0]
            else:
                tr.gene_id = attrib_dict['gene_id']
            
            tr_dict[tr_id] = tr
            tr_chrom_dict[chrom_name].append(tr)
            
        elif feature_type == 'exon':
            # create new exon
            ex = genome.coord.Coord(chrom_dict[chrom_name],
                                    start, end, strand)
            
            if strip_version:
                ex.gene_id = attrib_dict['gene_id'].split(".")[0]
            else:
                ex.gene_id = attrib_dict['gene_id']

            tr_id = attrib_dict['transcript_id']
            if strip_version:
                tr_id = tr_id.split(".")[0]
            ex.transcript_id = tr_id

            # record transcript => exons mapping
            if tr_id in tr_ex_dict:
                tr_ex_dict[tr_id].append(ex)
            else:
                tr_ex_dict[tr_id] = [ex]

        elif feature_type == 'gene':
            gene = genome.gene.Gene()
            gene.chrom = chrom_dict[chrom_name]
            gene.start = start
            gene.end = end
            gene.strand = strand
            
            gene.gene_id = attrib_dict['gene_id']
            if strip_version:
                gene.gene_id = gene.gene_id.split(".")[0]

            gene.gene_name = attrib_dict['gene_name']
                        
            gene_dict[gene.gene_id] = gene
            gene_chrom_dict[chrom_name].append(gene)


    # now that all exons and genes have been parsed,
    # go through and add exons to transcripts, and transcripts
    # to genes
    for tr in tr_dict.values():
        if tr.gene_id in gene_dict:
            gene_dict[tr.gene_id].add_transcript(tr)
        else:
            # some genes (pseudogenes?) are not in GFF file
            gene = genome.gene.Gene()
            gene.gene_id = tr.gene_id
            gene_dict[gene_id] = gene
            gene.add_transcript(tr)

        if tr.transcript_id in tr_ex_dict:
            for exon in tr_ex_dict[tr.transcript_id]:
                tr.add_exon(exon)
        else:
            raise ValueError("file does not contain any exons "
                             "for transcript %s\n" % 
                             tr.transcript_id)
            
    f.close()

    # sort lists of genes and transcripts in chrom dict
    for chrom_name in gene_chrom_dict.keys():
        genome.coord.sort_coords(gene_chrom_dict[chrom_name])
        genome.coord.sort_coords(tr_chrom_dict[chrom_name])
    
    return gene_dict, tr_dict, gene_chrom_dict, tr_chrom_dict



