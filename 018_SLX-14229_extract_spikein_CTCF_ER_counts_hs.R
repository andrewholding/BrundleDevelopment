# TODO: Add comment
# 
# Author: giorgi01 & holdin01
###############################################################################


### Essentially the idea is to use spike in reads to normalize the Human reads
library(DiffBind)

setwd("/Volumes/FlyPeaks/flypeaks")
filename<-"Rdata/018_SLX-14229_dba_human_CTCF.rda"
if(!file.exists(filename)){
  dba<-dba(sampleSheet = "samplesheet/samplesheet_SLX14229_hs_CTCF.csv")
  save(dba,file=filename)
} else {
  load(filename)
}


### Ash's magic
op<-dba.overlap(dba, mode=DBA_OLAP_RATE)
png("plots/018_SLX-14229_diffbind_overlap_human_CTCF.png",w=1000,h=1000,p=30)
plot(op,type="o",lwd=3,xlab="Nr. Samples",main="Binding Overlaps",ylab="Nr. Peaks",pch=16)
grid()
dev.off()



filename<-"Rdata/018_SLX-14229_dba.counts_human_CTCF.rda"
if(!file.exists(filename)){
  dba <- dba.count(dba, minOverlap = 5, summits=200)
  
  ### If peaks==NULL the ‘score’, ‘filter’, and ‘summits’ parameters are honored, updating the global binding matrix
  ### without re-counting in the cases of ‘score’ and ‘filter’, and only counting after re-centering
  ### in the case of "summits"
  ### Asking for DBA_SCORE_READS provides raw read counts from ChIP only (the simplest score)
  dba <- dba.count(dba, peaks=NULL, score=DBA_SCORE_READS)
  
  png("plots/018_SLX-14229_dba.counts_CTCF.png")
  plot(dba)
  dev.off()
  
  save(dba,file=filename)
} else {
  load(filename)
}


### Extract the peakset as a matrix
hsconsensus <- dba.peakset(dba, bRetrieve = T, DataType = DBA_DATA_FRAME)
save(hsconsensus,file="Rdata/018_SLX-14229_hsconsensus_CTCF.rda")
saveRDS(hsconsensus,file="Rdata/018_SLX-14229_hsconsensus_CTCF.rds")


png("plots/018_SLX-14229_MAplot_human_CTCF.png")
dba.plotMA(dba.analyze(dba))
dev.off()

rm(list=ls())

#Do the same again for ER peaks


setwd("/Volumes/FlyPeaks/flypeaks")
filename<-"Rdata/018_SLX-14229_dba_human_ER.rda"
if(!file.exists(filename)){
  dba<-dba(sampleSheet = "samplesheet/samplesheet_SLX14229_hs_ER.csv")
  save(dba,file=filename)
} else {
  load(filename)
}


### Ash's magic
op<-dba.overlap(dba, mode=DBA_OLAP_RATE)
png("plots/018_SLX-14229_diffbind_overlap_human_ER.png",w=1000,h=1000,p=30)
plot(op,type="o",lwd=3,xlab="Nr. Samples",main="Binding Overlaps",ylab="Nr. Peaks",pch=16)
grid()
dev.off()



filename<-"Rdata/018_SLX-14229_dba.counts_human_ER.rda"
if(!file.exists(filename)){
  dba <- dba.count(dba, summits=200)
  
  ### If peaks==NULL the ‘score’, ‘filter’, and ‘summits’ parameters are honored, updating the global binding matrix
  ### without re-counting in the cases of ‘score’ and ‘filter’, and only counting after re-centering
  ### in the case of "summits"
  ### Asking for DBA_SCORE_READS provides raw read counts from ChIP only (the simplest score)
  dba <- dba.count(dba, peaks=NULL, score=DBA_SCORE_READS)
  
  png("plots/018_SLX-14229_dba.counts_human_ER.png")
  plot(dba)
  dev.off()
  
  save(dba,file=filename)
} else {
  load(filename)
}


### Extract the peakset as a matrix
hsconsensus <- dba.peakset(dba, bRetrieve = T, DataType = DBA_DATA_FRAME)
save(hsconsensus,file="Rdata/018_SLX-14229_hsconsensus_ER.rda")
saveRDS(hsconsensus,file="Rdata/018_SLX-14229_hsconsensus_ER.rds")


png("plots/018_SLX-14229_MAplot_human_ER.png")
dba.plotMA(dba.analyze(dba))
dev.off()
