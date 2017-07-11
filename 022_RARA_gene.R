source("https://bioconductor.org/biocLite.R")
biocLite("Gviz")

library(Gviz)



bamFile <- "./ENCODE/ENCFF591QUF.bam.fq.hg19.bam"
dTrack3 <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "ENCFF783NXW", window = -1, chromosome = "chr17", col="green", )
plotTracks(dTrack3, chromosome = "chr17", from =  38470000, to = 38520000)


bamFile <- "./ENCODE/ENCFF783NXW.bam.fq.hg19.bam"
dTrack4 <- DataTrack(range = bamFile, genome = "hg19", type = "l", name = "ENCFF783NXW", window = -1, chromosome = "chr17", col="green", )
plotTracks(dTrack4, chromosome = "chr17", from =  38470000, to = 38520000)

ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chr17")
plotTracks(ideoTrack, from = 38470000, to = 38520000)

axisTrack <- GenomeAxisTrack()
plotTracks(axisTrack, from = 38470000, to = 38520000,  chromosome = "chr17", scale = 0.5)

#Generate annotation - https://support.bioconductor.org/p/69602/
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
txTranscripts_v1 <- GeneRegionTrack(txdb_hg19, genome="hg19", chromosome="chr17", showId=TRUE, geneSymbol=TRUE, name="UCSC")

library(org.Hs.eg.db)
symbols <- unlist(mapIds(org.Hs.eg.db, gene(txTranscripts_v1), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(txTranscripts_v1) <- symbols[gene(txTranscripts_v1)]

plotTracks(txTranscripts_v1, from = 38470000, to = 38520000)



plotTracks(list(ideoTrack, axisTrack,dTrack3, dTrack4, txTranscripts_v1), from =  38470000, to = 38520000)
