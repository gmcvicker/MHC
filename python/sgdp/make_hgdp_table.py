import sys


# file containing list of HGDP samples in SGDP project
SGDP_HGDP_TABLE = "/home/gm114/data/SGDP/samples/FullyPublicHGDP/Summary/Summary_info.txt"
SGDP_GENERAL_TABLE = "/home/gm114/data/SGDP/samples/FullyPublicGeneral/Summary/Summary_info.txt"
HGDP_TABLE = "/home/gm114/data/hgdp/ceph-hgdp/hgdp-ceph-unrelated.out"


# read HGDP identifiers that are in SGDP
f = open(SGDP_HGDP_TABLE)
header = f.readline()

sgdp_panel = {}
cteam_id = {}
illumina_id = {}

for line in f:
    words = line.split("\t")
    hgdp_id = words[3]
    panel = words[1]

    sgdp_panel[hgdp_id] = panel
    illumina_id[hgdp_id] = words[2]
    cteam_id[hgdp_id] = words[5]

f.close()


# read other SGDP samples (includes many non-HGDP)
f = open(SGDP_GENERAL_TABLE)
header = f.readline()

for line in f:
    words = line.split("\t")
    sample_id = words[3]
    panel = words[1]

    if sample_id not in sgdp_panel:
        sgdp_panel[sample_id] = panel
        illumina_id[sample_id] = words[2]
        cteam_id[sample_id] = words[5]
    
    


f = open(HGDP_TABLE)

sys.stdout.write("\t".join(["HGDP.ID", "GENDER", "POPULATION", "GEOGRAPHIC.ORIGIN", 
                            "REGION", "SGDP.PANEL", "ILLUMINA.ID", "C.TEAM.ID"]) + "\n")

header = f.readline()

for line in f:
    words = line.rstrip().split("\t")
    hgdp_id = words[0]
    if hgdp_id in sgdp_panel:
        words.extend([sgdp_panel[hgdp_id], illumina_id[hgdp_id], cteam_id[hgdp_id]])
    else:
        words.extend([".", ".", "."])
        
    sys.stdout.write("\t".join(words) + "\n")






    
