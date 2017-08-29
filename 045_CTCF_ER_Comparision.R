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

#Sequence
#library(BSgenome.Hsapiens.UCSC.hg19)
#seqTrack <- SequenceTrack(Hsapiens)
#plotTracks(seqTrack, from = windowStart, to = windowEnd-950,  chromosome = chromosome)

#Motif
#Full Genome Output from Homer
#AAGGTCAC = ER
#CRGCAGAGGGCR = CTCF
#    chr9    97545085        97545094        8-GTAATTACTT    9       -
#    chr9    97545192        97545199        1-AAGGTCAC      6       +
#    chr9    97545279        97545286        1-AAGGTCAC      5       -
#    chr9    97545386        97545395        4-GGGTGACCTT    8       +
#    chr9    97545388        97545395        1-AAGGTCAC      5       -
#    chr9    97545401        97545408        1-AAGGTCAC      7       +
#    chr9    97545406        97545415        9-CAGAGTAGCC    9       +
#    chr9    97545504        97545513        1-CACCAGAGGG    8       +
#    chr9    97545504        97545515        1-CRGCAGAGGGCR  12      +
#    chr9    97545509        97545516        3-GAGGGCGC      9       +
#    chr9    97545538        97545545        1-AAGGTCAC      5       -
#    chr9    97545547        97545554        1-AAGGTCAC      7       +
#    chr9    97545547        97545556        4-GGGTGACCTT    8       -
#    chr9    97545843        97545852        7-TCCTGRCASA    12      +
#CTC site at 97545000+504
#ERE site at 97545000+192
#ERE site at 97545000+279
#ERE site at 97545000+388
#ERE site at 97545000+401
#ERE site at 97545000+538
#ERE site at 97545000+547
motifTrack <- AnnotationTrack(start = c(windowStart+504, windowStart+192,
                                        windowStart+279, windowStart+388,
                                        windowStart+401, windowStart+538,
                                        windowStart+547),
                           width = c(12,8,8,8,8,8,8), chromosome = "chr9", 
                           strand = "*", 
                           id = c("CTCF","ERE","ERE","ERE","ERE","ERE","ERE"),
                           genome = "hg19", name = "")
feature(motifTrack)<-c("CTCF","ERE","ERE","ERE","ERE","ERE","ERE")
plotTracks(motifTrack,featureAnnotation = "id", ERE = "blue", CTCF = "purple")



pdf("plots/045_RARA_locus_min.pdf",pointsize=3,width=8,height=2)
plotTracks(list(ideoTrack,dTrack1b,dTrack2a,CTCFb,motifTrack),
           from =  windowStart, to = windowEnd, #featureAnnotation = "id"
            ERE = "blue", CTCF = "purple")
dev.off()



