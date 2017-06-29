library(DiffBind)

setwd("/Volumes/FlyPeaks/flyPeaks")

dir.create("./plots")
dir.create("./Rdata")

### Mouse Peaks

dba<-dba(sampleSheet = "samplesheet/samplesheet_SLX12998.csv")

op<-dba.overlap(dba, mode=DBA_OLAP_RATE)

png("plots/006_SLX-12998_binding_overlaps.png",w=1000,h=1000,p=30)
plot(op,type="o",lwd=3,xlab="Nr. Samples",main="Binding Overlaps",ylab="Nr. Peaks",pch=16)
dev.off()

dba <- dba.count(dba,  summits=200)

png("plots/006_SLX-12998_heatmap.png",w=1000,h=1000,p=30)
plot(dba)
dev.off()

png("plots/006_SLX-12998_PCA.png",w=1000,h=1000,p=30)
dba.plotPCA(dba)
dev.off()

dba.contrast(dba)
dba_analyze<-dba.analyze(dba)

png("plots/006_SLX-12998_MA.png",w=1000,h=1000,p=30)
dba.plotMA(dba_analyze, contrast=1)
dev.off()


### Fly Peaks

dba_fly<-dba(sampleSheet = "samplesheet/samplesheet_SLX8047.csv")

op<-dba.overlap(dba_fly, mode=DBA_OLAP_RATE)

png("plots/006_SLX-8047_Binding_Overlaps.png",w=1000,h=1000,p=30)
plot(op,type="o",lwd=3,xlab="Nr. Samples",main="Binding Overlaps",ylab="Nr. Peaks",pch=16)
dev.off()

dba_fly <- dba.count(dba_fly,  summits=200)
png("plots/006_SLX-8047_heatmap.png",w=1000,h=1000,p=30)
plot(dba_fly)
dev.off()

png("plots/006_SLX-8047_PCA.png",w=1000,h=1000,p=30)
dba.plotPCA(dba_fly)
dev.off()

dba.contrast(dba_fly)
dba_fly_analyze<-dba.analyze(dba_fly)

png("plots/006_SLX-8047_MA.png",w=1000,h=1000,p=30)
dba.plotMA(dba_fly_analyze, contrast=1)
dev.off()

### CTCFy Peaks

dba_ctcf<-dba(sampleSheet = "samplesheet/samplesheet_SLX14229.csv")

op<-dba.overlap(dba_ctcf, mode=DBA_OLAP_RATE)

png("plots/006_SLX-14229_Binding_Overlaps.png",w=1000,h=1000,p=30)
plot(op,type="o",lwd=3,xlab="Nr. Samples",main="Binding Overlaps",ylab="Nr. Peaks",pch=16)
dev.off()

dba_ctcf <- dba.count(dba_ctcf,  summits=200) #Less reps so less overlap
png("plots/006_SLX-14229_heatmap.png",w=1000,h=1000,p=30)
plot(dba_ctcf)
dev.off()

png("plots/006_SLX-14229_PCA.png",w=1000,h=1000,p=30)
dba.plotPCA(dba_ctcf)
dev.off()

dba.contrast(dba_ctcf)
dba_ctcf_analyze<-dba.analyze(dba_ctcf)

png("plots/006_SLX-14229_MA.png",w=1000,h=1000,p=30)
dba.plotMA(dba_ctcf_analyze, contrast=1)
dev.off()

