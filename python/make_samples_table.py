import sys


BASE_DIR = "/home/gm114/data/SGDP/samples"

# Complete list of C-team samples
CTEAM_LIST = "%s/hetfa.dblist" % BASE_DIR

EMBARGO_LEVEL = ["FullyPublicHGDP",
                 "FullyPublicGeneral",
                 "SignedLetterDiRienzo",
                 "SignedLetterGeneral",
                 "SignedLetterTishkoff",
                 "SignedLetterIndia"]


def read_cteam_index():
    f = open(CTEAM_LIST)
    i = 0
    
    cteam_idx_dict = {}
    
    for line in f:
        words = line.rstrip().split()

        cteam_id = words[0]

        cteam_idx_dict[cteam_id] = i

        i += 1

    f.close()

    return cteam_idx_dict
        


    
def main():
    cteam_idx_dict = read_cteam_index()

    first = True
    
    for embargo_level in EMBARGO_LEVEL:
        filename = "%s/%s/Summary/Summary_info.txt" % (BASE_DIR, embargo_level)

        sys.stderr.write("%s\n" % filename)

        f = open(filename)

        header = f.readline()
        
        if first:
            # only write header for first file
            # add CTeam_Index as column
            header = header.rstrip() + "\tCTeam_Index"
            
            sys.stdout.write(header + "\n")
            first = False

        for line in f:
            line = line.rstrip()
            words = line.split("\t")

            cteam_id = words[5]
            if cteam_id in cteam_idx_dict:
                cteam_idx = cteam_idx_dict[cteam_id]
            else:
                cteam_idx = "NA"

            sys.stdout.write("%s\t%d\n" % (line, cteam_idx))

        f.close()


    
main()
