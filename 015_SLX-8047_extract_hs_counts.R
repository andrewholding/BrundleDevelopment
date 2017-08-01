# Extract counts
# 
# Author: giorgi01 & holdin01
###############################################################################



library(DiffBind)

setwd("/Volumes/FlyPeaks/flypeaks")

### Number of Drosophila reads in peaks
load("Rdata/014_SLX-8047_dmconsensus.rda")

### Set-up the human DBA object
filename<-"Rdata/015_SLX-8047_dba_human.rda"
if(!file.exists(filename)){
	dba<-dba(sampleSheet = "samplesheet/samplesheet_SLX8047_hs.csv")
	save(dba,file=filename)
} else {
	load(filename)
}

### Correlation between samples
png("plots/015_SLX-8047_diffbind_samplesheet_human.png",w=1000,h=1000,p=30)
plot(dba)
dev.off()

### Compute binding overlaps and co-occupancy statistics
dbaoverlap<-dba.overlap(dba, mode=DBA_OLAP_RATE)
png("plots/015_SLX-8047_diffbind_overlap_human.png",w=1000,h=1000,p=30)
plot(dbaoverlap,type="o",lwd=3,xlab="Nr. Samples",main="Binding Overlaps",ylab="Nr. Peaks",pch=16)
grid()
dev.off()


### Count reads in binding sites intervals
filename<-"Rdata/015_SLX-8047_dbacounts_human.rda"
if(!file.exists(filename)){
	dbacounts <- dba.count(dba, summits = 200)
	save(dbacounts,file=filename)
} else {
	load(filename)
}
png("plots/015_SLX-8047_diffbind_dbacounts_human.png",w=1000,h=1000,p=30)
plot(dbacounts)
dev.off()
png("plots/015_SLX-8047_diffbind_dbacounts_nocorr_human.png",w=1000,h=1000,p=30)
plot(dbacounts,correlations=F,maxval=6)
dev.off()



### Perform differential binding affinity analysis
filename<-"Rdata/015_SLX-8047_dbanalysis.rda"
if(!file.exists(filename)){
	dbanalysis <- dba.analyze(dbacounts)
	save(dbanalysis,file=filename)
} else {
	load(filename)
}

##########################################################
### Ok now we prove the concept with a simple T-test
dbacounts2 <- dba.count(dbacounts, peaks=NULL, score=DBA_SCORE_READS)
hsconsensus<-dba.peakset(dbacounts2, bRetrieve = T, DataType = DBA_DATA_FRAME)
head(hsconsensus) # contains the reads in peaks data
head(dmconsensus) # reads in peaks for drosophila melanogaster
# Input has been used already as a subtraction sample for reads
saveRDS(hsconsensus,file="Rdata/015_SLX-8047_hsconsensus.rds")
saveRDS(dmconsensus,file="Rdata/015_SLX-8047_dmconsensus.rds")


### Plot MA

dba.contrast(dbanalysis)
png("plots/015_SLX-8047_plotMA_defaults.png")
dba.plotMA(dbanalysis,contrast=1)
dev.off()

### Normalize by reads in peakset
hscounts<-hsconsensus[,-c(1:3)]
dmcounts<-dmconsensus[,-c(1:3)]
hssums<-apply(hscounts,2,sum)/1e6
dmsums<-apply(dmcounts,2,sum)/1e6
names(dmsums)<-names(hssums)<-gsub("X","",names(hssums))

png("plots/015_SLX-8047_readsInPeaks.png",w=1000,h=1000,p=40)
plot(hssums*1e3,dmsums*1e3,pch=16,xlab="Human (K)",ylab="Drosophila (K)",col="lightgrey",main="Total Reads in Peaks",xlim=c(0,max(hssums*1e3)))
grid()
text(hssums*1e3,dmsums*1e3,labels=names(hssums))
mtext("a=Fulvestrant, b=None",cex=0.8)
dev.off()



hsnorm<-as.matrix(hscounts)


# Multiply by the drosophila normalization factor
#hsnorm2<-hsnorm
#dmnormfactor<-dmsums/mean(dmsums)
#for(i in 1:ncol(hsnorm)){
#  hsnorm2[,i]<-hsnorm[,i]*dmnormfactor[i]
#}

# Divide by the drosophila normalization factor
hsnorm3<-hsnorm
dmnormfactor<-dmsums/mean(dmsums)
for(i in 1:ncol(hsnorm)){
  hsnorm3[,i]<-hsnorm[,i]/dmnormfactor[i]
}

png("plots/015_SLX-8047_readsInPeaks_dmDivide.png",w=1000,h=1000,p=40)
hssums3<-apply(hsnorm3,2,sum)/1e6
names(hssums3)<-gsub("X","",names(hssums))
plot(hssums3*1e3,dmsums*1e3,pch=16,xlab="Human (K / (dm reads in peaks))",ylab="Drosophila (K)",col="lightgrey",main="Total Reads in Peaks",xlim=c(0,max(hssums*1e3)))
grid()
text(hssums3*1e3,dmsums*1e3,labels=names(hssums3))
mtext("a=Fulvestrant, b=None",cex=0.8)
dev.off()












