#source("https://bioconductor.org/biocLite.R")
#biocLite("Gviz")

library(Gviz)
windowStart<-38470000
windowEnd <- 38520000
chromosome<-"chr17"


#windowStart<-43780391 -1000
#windowEnd<-43788644 + 1000
#chromosome<-"chr21"

setwd("/Volumes/FlyPeaks/FlyPeaks")

#-ve in sample sheet for SLX-14229 is ICI, +ve is without.

bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D701_D503.HJJL7BBXX.s_8.r_1.fq.gz.bam"
dTrack1a <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "ICI1", window = -1, chromosome = chromosome, col="red", ylim=c(0,225),showAxis=F)
plotTracks(dTrack1a, chromosome = chromosome, from =  windowStart, to = windowEnd)


bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D702_D504.HJJL7BBXX.s_8.r_1.fq.gz.bam"
dTrack2a <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "ICI2", window = -1, chromosome = chromosome, col="red", ylim=c(0,225),showAxis=F)
plotTracks(dTrack2a, chromosome = chromosome, from =  windowStart, to = windowEnd)


bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D703_D502.HJJL7BBXX.s_8.r_1.fq.gz.bam"
dTrack3a <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "ICI3", window = -1, chromosome = chromosome, col="red",ylim=c(0,225),showAxis=F )
plotTracks(dTrack3a, chromosome = chromosome, from =  windowStart, to = windowEnd)





bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D701_D501.HJJL7BBXX.s_8.r_1.fq.gz.bam"
dTrack1b <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "Ctrl1", window = -1, chromosome = chromosome, col="blue",ylim=c(0,225),showAxis=F )
plotTracks(dTrack1b, chromosome = chromosome, from =  windowStart, to = windowEnd)


bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D702_D502.HJJL7BBXX.s_8.r_1.fq.gz.bam"
dTrack2b <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "Ctrl2", window = -1, chromosome = chromosome, col="blue", ylim=c(0,225),showAxis=F)
plotTracks(dTrack2b, chromosome = chromosome, from =  windowStart, to = windowEnd)


bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D703_D504.HJJL7BBXX.s_8.r_1.fq.gz.bam"
dTrack3b <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "Ctrl3", window = -1, chromosome = chromosome, col="blue",ylim=c(0,225),showAxis=F )
plotTracks(dTrack3b, chromosome = chromosome, from =  windowStart, to = windowEnd)

#CTCF with ICI

bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D701_D504.HJJL7BBXX.s_8.r_1.fq.gz.bam"
CTCFa <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "CTCF", window = -1, chromosome = chromosome, col="purple",ylim=c(0,225),showAxis=F)
plotTracks(CTCFa, chromosome = chromosome, from =  windowStart, to = windowEnd)

#CTCF without ICI

bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D702_D503.HJJL7BBXX.s_8.r_1.fq.gz.bam"
CTCFb <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "CTCF", window = -1, chromosome = chromosome, col="purple",showAxis=F )
plotTracks(CTCFb, chromosome = chromosome, from =  windowStart, to = windowEnd)


bamFile <- "./ENCODE/ENCFF591QUF.bam.fq.hg19.bam"
Encode1 <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "ENCODE", window = -1, chromosome = chromosome, col="green",showAxis=F)


bamFile <- "./ENCODE/ENCFF783NXW.bam.fq.hg19.bam"
Encode2 <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "ENC2", window = -1, chromosome = chromosome, col="green", showAxis=F)
plotTracks(Encode2, chromosome = chromosome, from =  windowStart, to = windowEnd)

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


png("plots/022_RARA_locus_full.png",w=1000,h=1000,p=9)
plotTracks(list(ideoTrack, axisTrack,dTrack1a,dTrack2a,dTrack3a,dTrack1b,dTrack2b,dTrack3b,CTCFa,CTCFb,Encode1, Encode2, txTranscripts_v1), from =  windowStart, to = windowEnd)
dev.off()

png("plots/022_RARA_locus_clear.png",w=1000,h=1000,p=9)
plotTracks(list(ideoTrack, axisTrack,dTrack1a,dTrack2a,dTrack3a,dTrack1b,dTrack2b,dTrack3b,CTCFa,Encode1,  txTranscripts_v1), from =  windowStart, to = windowEnd)
dev.off()


pdf("plots/022_RARA_locus_min.pdf",pointsize=3)
plotTracks(list(ideoTrack, axisTrack,dTrack1a,dTrack1b,CTCFa,  txTranscripts_v1), from =  windowStart, to = windowEnd)
dev.off()


dTrack1a_scaled <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "ICI1", window = -1, chromosome = chromosome, col="red", ylim=c(0,96),showAxis=F)
CTCFa_scaled <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "CTCF", window = -1, chromosome = chromosome, col="purple",ylim=c(0,96),showAxis=F)
pdf("plots/022_RARA_locus_min_scaled.pdf",pointsize=3)
plotTracks(list(ideoTrack, axisTrack,dTrack1a_scaled,dTrack1b,CTCFa_scaled,  txTranscripts_v1), from =  windowStart, to = windowEnd)
dev.off()


