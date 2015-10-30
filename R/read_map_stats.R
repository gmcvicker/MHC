library(RColorBrewer)

samples.tab <- read.table("~/data/MHC/Pilot_Samples_MiSeq_151005.tsv", header=T, sep="\t")
tab <- read.table("~/data/MHC/MiSeq_151005/bwa/read_map_stats.txt.gz", header=T)

pal <- brewer.pal(nlevels(small.tab$MHC.STATUS), "Blues")
pal[levels(small.tab$MHC.STATUS) == "NON-MHC"] <- "red"

loci <- c("DPB1", "DQA1", "DQB1", "DRB", "HLA-A", "HLA-B", "HLA-C")
short.ids <- c("449", "476", "915", "928", "932", "H991", 
               "P1038", "P107", "P1103", "P253", "P496", "P652")

for(locus in loci) {
  pdf(paste("read_map_counts.", locus, ".pdf", sep=""),
      width=15, height=10)
  
  par(mfrow=c(3, length(short.ids)/3))
  
  for(short.id in short.ids) {
      species <- samples.tab[(short.id == samples.tab$SHORT.ID), "SHORT.SPECIES"]
      
      small.tab <- tab[(tab$LOCUS == locus) & (tab$SHORT.ID==short.id),]
      
      clrs <- pal[small.tab$MHC.STATUS]
      barplot(small.tab$N.MAPPED, col=clrs, border="black", 
              xlab="num mapped reads", las=1,
              ylim=c(0, 20), cex.axis=0.5,
              horiz=T, names.arg=substr(small.tab$SEQ.NAME, 0, 8),
              main=paste(locus, " - ", short.id, " (", species, ")", sep=""))
      
      legend("topright", as.character(levels(small.tab$MHC.STATUS)), fill=pal)
  }
  
  dev.off()
}


