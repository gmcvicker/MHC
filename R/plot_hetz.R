
#tab <- read.table("/home/gm114/data/SGDP/mhc_hetz.txt.gz",
#                  header=F)


## classic.loci.tab <- read.table("/home/gm114/data/GENCODE/hg19/mhc_classical_genes.txt",
##                                header=F, col.names=c("GENE.ID", "GENE.NAME", "CHROM",
##                                              "GENE.NUM", "START", "END", "STRAND"))


hla.colors <- data.frame(TYPE=c("Class II Classical",
                             "Class II Non-classical",
                             "Class I Non-classical",
                             "Class I Classical"),
                         COLOR=c("green", "blue", "red", "pink"))

classic.loci.tab <- read.table("/home/gm114/data/GENCODE/hg19/mhc_gene_types.txt",
                               sep="\t", header=F,
                               col.names=c("GENE.ID", "GENE.NAME", "CHROM",
                                   "GENE.NUM", "START", "END", "STRAND",
                                   "TYPE"))
classic.loci.tab <- merge(classic.loci.tab, hla.colors, by="TYPE", all.x=TRUE)


regions <- c("WestEurasia", "SouthAsia", "EastAsia", "CentralAsiaSiberia",
             "America", "Africa", "Oceania")



gene.tab <- read.table("/home/gm114/data/GENCODE/hg19/mhc_all_genes.txt", header=F,
                       col.names=c("GENE.ID", "GENE.NAME", "CHROM", "GENE.NUM",
                           "START", "END", "STRAND"))

exon.tab <- read.table("/home/gm114/data/GENCODE/hg19/mhc_all_exons.txt", header=F,
                       col.names=c("GENE.ID", "GENE.NAME", "CHROM", "EXON.NUM",
                           "START", "END", "STRAND"))


# add HLA type labels to genes and exons
gene.tab <- merge(gene.tab, classic.loci.tab[,c("GENE.ID", "TYPE", "COLOR")],
                   by="GENE.ID", all.x=T)

exon.tab <- merge(exon.tab, classic.loci.tab[,c("GENE.ID", "TYPE", "COLOR")],
                   by="GENE.ID", all.x=T)



plot.genes <- function(gene.tab, exon.tab,
                       y.mid=-0.01, exon.height=0.005) {
    
    y.mid.fwd <- y.mid + exon.height
    y.mid.rev <- y.mid - exon.height
    
    # draw introns first

    y.mid <- rep(y.mid, nrow(gene.tab))
    y.mid[gene.tab$STRAND == 1] <- y.mid.fwd
    y.mid[gene.tab$STRAND == -1] <- y.mid.rev    

    color <- as.character(gene.tab$COLOR)
    color[is.na(color)] <- "grey50"
    
    segments(x0=gene.tab$START, x1=gene.tab$END,
             y0=y.mid, y1=y.mid,
             col=color)
    
    y.mid <- rep(y.mid, nrow(exon.tab))
    y.mid[exon.tab$STRAND == 1] <- y.mid.fwd
    y.mid[exon.tab$STRAND == -1] <- y.mid.rev    
    
    color <- as.character(exon.tab$COLOR)
    color[is.na(color)] <- "grey50"
    
    # draw exons
    rect(xleft=exon.tab$START, ybottom=y.mid-exon.height/2,
         xright=exon.tab$END, ytop=y.mid+exon.height/2,
         col=color, border=color)
    
}



## # build table of heterozygosity for each of the populations
## hetz.tab <- NULL

## for(region in regions) {
##     cat(region, "\n")
##     filename <- paste("/home/gm114/data/SGDP/mhc_hetz.", region, ".txt.gz", sep="")
##     tab <- read.table(filename)

##     if(is.null(hetz.tab)) {
##         hetz.tab <- data.frame(POS=tab$V1)
##     }

##     hetz.tab[[region]] <- tab$V2
## }


hetz.tab <- NULL

species.tab <- read.table("~/data/great_ape_genome_project/species.txt")
species <- as.character(species.tab$V1)

for(sp in species) {
    cat(sp, "\n")
    filename <- paste("~/data/great_ape_genome_project/mhc_hetz.", sp, ".txt.gz",
                      sep="")
    tab <- read.table(filename)

    if(is.null(hetz.tab)) {
        hetz.tab <- data.frame(POS=tab$V1)
    }
    hetz.tab[[sp]] <- tab$V2
}





mhc.start <- 28000000
pos <- hetz.tab$POS + mhc.start - 1
mhc.end <- max(pos)


types <- names(hetz.tab)[2:ncol(hetz.tab)]

win.size <- 1000
smooth.hetz.tab <- data.frame(hetz.tab)
for(type in types) {
    cat(type, "\n")
    
    smooth.hetz.tab[,type] <- filter(hetz.tab[,type],
                                       rep(1, win.size)) / win.size
}


# write smoothed table, since takes some time to create
#write.table(smooth.hetz.tab, quote=F, sep="\t", row.names=F, col.names=T,
#    file=paste("/home/gm114/data/SGDP/mhc_hetz.combined.smooth", win.size, ".txt", sep=""))

# smooth.hetz.tab <- read.table(paste("/home/gm114/data/SGDP/mhc_hetz.combined.smooth", win.size, ".txt.gz", sep=""), header=T)


segments.tab <-
    data.frame(NAME=c("wholeregion", "1", "2", "3",
               as.character(classic.loci.tab$GENE.NAME)),
        START=c(mhc.start, 29000000, 31000000, 32000000,
                   classic.loci.tab$START-1000),
               END=c(mhc.end,   30100000, 31500000, 33100000,
                   classic.loci.tab$END+1000))




library(RColorBrewer)
type.colors <- brewer.pal(length(types), "Set1")

for(seg in 1:nrow(segments.tab)) {
    segname <- segments.tab$NAME[seg]
    
    # png(paste("mhc_hetz.", segname, ".png", sep=""),  width=1000, height=500)
    pdf(paste("mhc_hetz.", segname, ".pdf", sep=""),  width=10, height=5)

    xlim <- c(segments.tab$START[seg], segments.tab$END[seg])
            
    plot(c(0), c(0), xlab="chr 6 position", type="n",
         xlim=xlim, ylim=c(-0.02, 0.06),
         ylab="nucleotide diversity (pi)")

    lines(x=c(mhc.start, mhc.end), y=c(0,0), col="black")

    
    plot.genes(gene.tab, exon.tab)
    
    mid <- (classic.loci.tab$END + classic.loci.tab$START)/2

    y <- rep(c(-0.005, -0.01, -0.015), nrow(classic.loci.tab))[1:nrow(classic.loci.tab)]
    
    text(mid, y=y, labels=as.character(classic.loci.tab$GENE.NAME),
         pos=1, col=as.character(classic.loci.tab$COLOR))
    
    # plot diversity for each geographic type
    for(i in 1:length(types)) {
        type <- types[i]
        lines(pos, smooth.hetz.tab[[type]], col=type.colors[i])
    }


    legend("topleft", legend=types, col=type.colors, lty=1, bg="white")
    
    dev.off()    
}





