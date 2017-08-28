#source("https://bioconductor.org/biocLite.R")
#biocLite("Gviz")

library(Gviz)
windowStart<-97545000
windowEnd <- 97546000
chromosome<-"chr9"


setwd("/Volumes/FlyPeaks/FlyPeaks")

#-ve in sample sheet for SLX-14229 is ICI, +ve is without.

bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D702_D504.HJJL7BBXX.s_8.r_1.fq.gz.bam"
dTrack2a <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "ICI", window = -1, chromosome = chromosome, col="red", ylim=c(0,225),showAxis=F)
plotTracks(dTrack2a, chromosome = chromosome, from =  windowStart, to = windowEnd)
##Keeep CONTROL

bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D701_D501.HJJL7BBXX.s_8.r_1.fq.gz.bam"
dTrack1b <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "Ctrl", window = -1, chromosome = chromosome, col="blue",ylim=c(0,225),showAxis=F )
plotTracks(dTrack1b, chromosome = chromosome, from =  windowStart, to = windowEnd)

#KEEP CTCF without ICI

bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D702_D503.HJJL7BBXX.s_8.r_1.fq.gz.bam"
CTCFb <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "CTCF", window = -1, chromosome = chromosome, col="purple",showAxis=F )
plotTracks(CTCFb, chromosome = chromosome, from =  windowStart, to = windowEnd)

ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chromosome)
plotTracks(ideoTrack, from = windowStart, to = windowEnd)

axisTrack <- GenomeAxisTrack()
plotTracks(axisTrack, from = windowStart, to = windowEnd,  chromosome = chromosome, scale = 0.5)

#Generate annotation - https://support.bioconductor.org/p/69602/
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
txTranscripts_v1 <- GeneRegionTrack(txdb_hg19, genome="hg19", chromosome=chromosome, showId=TRUE, geneSymbol=TRUE, name="UCSC")

library(org.Hs.eg.db)
symbols <- unlist(mapIds(org.Hs.eg.db, gene(txTranscripts_v1), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(txTranscripts_v1) <- symbols[gene(txTranscripts_v1)]

plotTracks(txTranscripts_v1, from = windowStart, to = windowEnd)




pdf("plots/045_RARA_locus_min.pdf",pointsize=3,width=8,height=2)
plotTracks(list(ideoTrack,dTrack1b,dTrack2a,CTCFb), from =  windowStart, to = windowEnd)
dev.off()



