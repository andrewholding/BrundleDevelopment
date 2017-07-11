#source("https://bioconductor.org/biocLite.R")
#biocLite("Gviz")

library(Gviz)
windowStart<-38470000
windowEnd <-38520000
chromosome<-"chr17"
setwd("/Volumes/FlyPeaks/FlyPeaks")


bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D702_D502.HJJL7BBXX.s_8.r_1.fq.gz.bam"
dTrack1 <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "D702_D502", window = -1, chromosome = chromosome, col="blue", )
plotTracks(dTrack1, chromosome = chromosome, from =  windowStart, to = windowEnd)


bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D703_D504.HJJL7BBXX.s_8.r_1.fq.gz.bam"
dTrack2 <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "D703_D504", window = -1, chromosome = chromosome, col="blue", )
plotTracks(dTrack2, chromosome = chromosome, from =  windowStart, to = windowEnd)


bamFile <- "SLX-14229_mmhs/blacklist_filtered/human/SLX-14229.D703_D504.HJJL7BBXX.s_8.r_1.fq.gz.bam"
dTrack2 <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "D703_D504", window = -1, chromosome = chromosome, col="blue", )
plotTracks(dTrack2, chromosome = chromosome, from =  windowStart, to = windowEnd)



bamFile <- "./ENCODE/ENCFF591QUF.bam.fq.hg19.bam"
dTrack3 <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "ENCFF783NXW", window = -1, chromosome = chromosome, col="green", )
plotTracks(dTrack3, chromosome = chromosome, from =  windowStart, to = windowEnd)


bamFile <- "./ENCODE/ENCFF783NXW.bam.fq.hg19.bam"
dTrack4 <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "ENCFF783NXW", window = -1, chromosome = chromosome, col="green", )
plotTracks(dTrack4, chromosome = chromosome, from =  windowStart, to = windowEnd)

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



plotTracks(list(ideoTrack, axisTrack,dTrack1,dTrack2,dTrack3, dTrack4, txTranscripts_v1), from =  windowStart, to = windowEnd)
