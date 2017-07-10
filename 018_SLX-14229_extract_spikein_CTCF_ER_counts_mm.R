# TODO: Add comment
# 
# Author: giorgi01 & holdin01
###############################################################################


### Essentially the idea is to use spike in reads to normalize the Human reads
library(DiffBind)

setwd("/Volumes/FlyPeaks/flypeaks")
filename<-"Rdata/018_SLX-14229_dba_mouse_CTCF.rda"
if(!file.exists(filename)){
  dba<-dba(sampleSheet = "samplesheet/samplesheet_SLX14229_mm_CTCF.csv")
  save(dba,file=filename)
} else {
  load(filename)
}


### Ash's magic
op<-dba.overlap(dba, mode=DBA_OLAP_RATE)
png("plots/018_SLX-14229_diffbind_overlap_mouse_CTCF.png",w=1000,h=1000,p=30)
plot(op,type="o",lwd=3,xlab="Nr. Samples",main="Binding Overlaps",ylab="Nr. Peaks",pch=16)
grid()
dev.off()



filename<-"Rdata/018_SLX-14229_dba.counts_mouse_CTCF.rda"
if(!file.exists(filename)){
  dba <- dba.count(dba, minOverlap = 5, summits=200)
    #5 as we want this really consistent and it doesn't change
  
  ### If peaks==NULL the ‘score’, ‘filter’, and ‘summits’ parameters are honored, updating the global binding matrix
  ### without re-counting in the cases of ‘score’ and ‘filter’, and only counting after re-centering
  ### in the case of "summits"
  ### Asking for DBA_SCORE_READS provides raw read counts from ChIP only (the simplest score)
  dba <- dba.count(dba, peaks=NULL, score=DBA_SCORE_READS)
  
  png("plots/018_SLX-14229_dba.counts_mouse_CTCF.png")
  plot(dba)
  dev.off()
  
  save(dba,file=filename)
} else {
  load(filename)
}


### Extract the peakset as a matrix
mmconsensus <- dba.peakset(dba, bRetrieve = T, DataType = DBA_DATA_FRAME)
save(mmconsensus,file="Rdata/018_SLX-14229_mmconsensus_CTCF.rda")

png("plots/018_SLX-14229_MAplot_mouse_CTCF.png")
dba.plotMA(dba.analyze(dba))
dev.off()

rm(list=ls())

#Do the same again for ER peaks


setwd("/Volumes/FlyPeaks/flypeaks")
filename<-"Rdata/018_SLX-14229_dba_mouse_ER.rda"
if(!file.exists(filename)){
  dba<-dba(sampleSheet = "samplesheet/samplesheet_SLX14229_mm_ER.csv")
  save(dba,file=filename)
} else {
  load(filename)
}


### Ash's magic
op<-dba.overlap(dba, mode=DBA_OLAP_RATE)
png("plots/018_SLX-14229_diffbind_overlap_mouse_ER.png",w=1000,h=1000,p=30)
plot(op,type="o",lwd=3,xlab="Nr. Samples",main="Binding Overlaps",ylab="Nr. Peaks",pch=16)
grid()
dev.off()



filename<-"Rdata/018_SLX-14229_dba.counts_mouse_ER.rda"
if(!file.exists(filename)){
  dba <- dba.count(dba, summits=200)
  
  ### If peaks==NULL the ‘score’, ‘filter’, and ‘summits’ parameters are honored, updating the global binding matrix
  ### without re-counting in the cases of ‘score’ and ‘filter’, and only counting after re-centering
  ### in the case of "summits"
  ### Asking for DBA_SCORE_READS provides raw read counts from ChIP only (the simplest score)
  dba <- dba.count(dba, peaks=NULL, score=DBA_SCORE_READS)
  
  png("plots/018_SLX-14229_dba.counts_mouse_ER.png")
  plot(dba)
  dev.off()
  
  save(dba,file=filename)
} else {
  load(filename)
}


### Extract the peakset as a matrix
mmconsensus <- dba.peakset(dba, bRetrieve = T, DataType = DBA_DATA_FRAME)
save(mmconsensus,file="Rdata/018_SLX-14229_mmconsensus_ER.rda")

png("plots/018_SLX-14229_MAplot_mouse_ER.png")
dba.plotMA(dba.analyze(dba))
dev.off()
