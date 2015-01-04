
tab <- read.table("/home/gm114/proj/MHC/python/mhc_hetz_10.txt.gz",
                  header=F)


mhc.start <- 28000000

pos <- tab$V1 + mhc.start
hetz <- tab$V2


win.size <- 10000
smooth.hetz <- filter(hetz, rep(1, win.size)) / win.size


png("mhc_hetz.png", width=1000, height=500)

plot(pos, smooth.hetz, xlab="chr 6 position", ylab="heterozygosity")

dev.off()
