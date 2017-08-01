# Extract counts
# 
# Author: giorgi01 & holdin01
###############################################################################


### Essentially the idea is to use Drosophila reads to normalize the Human reads
library(DiffBind)

setwd("/Volumes/FlyPeaks/flypeaks")
filename<-"Rdata/014_SLX-8047_dba_drosophila.rda"
if(!file.exists(filename)){
  dba<-dba(sampleSheet = "samplesheet/samplesheet_SLX8047_dm.csv")
  save(dba,file=filename)
} else {
  load(filename)
}


### Ash's magic
op<-dba.overlap(dba, mode=DBA_OLAP_RATE)
png("plots/014_SLX-8047_diffbind_overlap_drosophila.png",w=1000,h=1000,p=30)
plot(op,type="o",lwd=3,xlab="Nr. Samples",main="Binding Overlaps",ylab="Nr. Peaks",pch=16)
grid()
dev.off()



filename<-"Rdata/014_SLX-8047_dba.counts.rda"
if(!file.exists(filename)){
  dba <- dba.count(dba, minOverlap = 5, summits=200)
  
  ### If peaks==NULL the ‘score’, ‘filter’, and ‘summits’ parameters are honored, updating the global binding matrix
  ### without re-counting in the cases of ‘score’ and ‘filter’, and only counting after re-centering
  ### in the case of "summits"
  ### Asking for DBA_SCORE_READS provides raw read counts from ChIP only (the simplest score)
  dba <- dba.count(dba, peaks=NULL, score=DBA_SCORE_READS)
  
  png("plots/014_SLX-8047_dba.counts.png")
  plot(dba)
  dev.off()
  
  save(dba,file=filename)
} else {
  load(filename)
}


### Extract the peakset as a matrix
dmconsensus <- dba.peakset(dba, bRetrieve = T, DataType = DBA_DATA_FRAME)
save(dmconsensus,file="Rdata/014_SLX-8047_dmconsensus.rda")

# Total reads in peaks (aka measure of efficiency of ChIPSeq)
#readsInPeaks<-apply(dmconsensus[,-(1:3)],2,sum)
#dba.normfac <- median(readsInPeaks)/readsInPeaks

#
#### Not used:
#dba.con.var <- apply(dba.con.counts, 1, var)
#dba.con.sd <- apply(dba.con.counts, 1, sd)
#dba.con.cov <- dba.con.sd/apply(dba.con.counts, 1, mean)
#
#cormat <- matrix(nc=8, nr=8)
#for(i in 1:7){
#	for(j in (i+1):8){
#		cormat[i,j] <- round(cor(dba.con.counts[,i], dba.con.counts[,j], method="spearman"), 3)
#	}}
#rownames(cormat) <- colnames(cormat) <- colnames(dba.con.counts)
