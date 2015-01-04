
Heng has assembled terabytes of contigs from the SGDP.

Heng has 254 contigs that he has built from the public samples... rest could be done as well...

Heng's contig files are located here:

/n/data1/hms/genetics/reich/1000Genomes/lh3/C-team-mag/mag/

Contigs are a consensus sequence so base calls are pretty good.. 
Phased up to read length (100bp) level.

- GRC38 has 8-12 haplotypes, but they are not complete - 
   do not span entire region
- Heng recommends building database using the GRC38 contigs, mapping to this. 


1. Map all contigs to reference haplotypes to pull out set that 
   are MHC contigs. In addition to MHC haplotypes also include the entire 
   region of MHC + some flanking.

2. Map the "exons" in the reference contigs to the set of MHC contigs that 
   I have created to find best-matching allele.

Try BWA-mem, could also BLAST.

