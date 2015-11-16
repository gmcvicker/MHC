
tab <- read.table("~/data/MHC/MiSeq_151005/velvet/blast_contigs_summary.txt.gz", header=T)

samples.tab <- read.table("~/data/MHC/Pilot_Samples_MiSeq_151005.tsv", header=T, sep="\t")

hit.type <- rep(NA, nrow(tab))
hit.type[tab$N.MHC.HIT > 0] <- "MHC"
hit.type[tab$N.NON.MHC.HIT > 0] <- "NON-MHC"
hit.type[(tab$N.MHC.HIT > 0) & (tab$N.NON.MHC.HIT > 0)] <- "BOTH"
hit.type[(tab$N.NON.MHC.HIT > 20)] <- "REPEAT"
tab["HIT.TYPE"] <- as.factor(hit.type)


loci <- unique(tab$LOCUS[tab$LOCUS != "Undetermined"])
short.ids <- unique(tab$SHORT.ID[tab$SHORT.ID != "."])

n.contig.types <- nlevels(tab$HIT.TYPE)
n.types <- length(loci) * length(short.ids) * n.contig.types
df <- data.frame(LOCUS=rep(NA, n.types), 
                 SHORT.ID=rep(NA, n.types),
                 CONTIG.TYPE=rep(NA, n.types),
                 N.CONTIGS=rep(NA, n.types))


i <- 0
for(locus in loci) {
  for(short.id in short.ids) {
    for(hit.type in levels(tab$HIT.TYPE)) {
        f <- (tab$LOCUS == locus) & (tab$SHORT.ID == short.id) & (tab$HIT.TYPE == hit.type)
        i <- i + 1
        df[i, "LOCUS"] <- locus
        df[i, "SHORT.ID"] <- short.id
        df[i, "CONTIG.TYPE"] <- hit.type
        df[i, "N.CONTIGS"] <- sum(f)
    }
  }
}


library(ggplot2)

pdf("assembled_contig_summary.pdf", width=8, height=8)

ggplot(df, aes(CONTIG.TYPE)) + 
  geom_bar(aes(CONTIG.TYPE, N.CONTIGS, fill=CONTIG.TYPE), stat="sum") +
  facet_grid(LOCUS ~ SHORT.ID) + ylab("number of contigs") + xlab("contig type") +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

dev.off()

pdf("assembled_contig_summary.thresh10.pdf", width=8, height=8)

df.threshold <- data.frame(df)
df.threshold$N.CONTIGS[df.threshold$N.CONTIGS > 10] <- 10
ggplot(df.threshold, aes(CONTIG.TYPE)) + 
  geom_bar(aes(CONTIG.TYPE, N.CONTIGS, fill=CONTIG.TYPE), stat="sum") +
  facet_grid(LOCUS ~ SHORT.ID) + ylab("number of contigs") + xlab("contig type") +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

dev.off()
