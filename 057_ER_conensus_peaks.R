
library(DiffBind)

setwd("/Volumes/FlyPeaks/FlyPeaks")

if(!file.exists("Rdata/057_dba.rda")) {
    dba<-dba(sampleSheet="samplesheet/samplesheet_SLX12998_hs.csv")
    dba<-dba.count(dba,  minOverlap=3, bParallel=FALSE)
    dbaPeakset<-dba.peakset(dba,  bRetrieve = T, DataType = DBA_DATA_FRAME)
    save(dba, file="Rdata/057_dba.rda")
} else {
    load(file="Rdata/057_dba.rda")
}

write.table(dbaPeakset[,1:3],"bed/ERConsensus.bed", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

setwd("./bed")
system("bedtools slop -i ERConsensus.bed -b 500 -g hg19.chrom.sizes > ERConsensus500.bed")
setwd("../SLX-14438_merged/peaks")
system("bedtools subtract -A -a CTCF_union.bed -b ../../bed/ERConsensus500.bed > CTCF_ERConsensusremoved.bed")
