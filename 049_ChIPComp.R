
#source("https://bioconductor.org/biocLite.R")
#biocLite()

biocLite("ChIPComp")
library(ChIPComp)


#Example

confs=makeConf(system.file("extdata", "conf.csv", package="ChIPComp"))
conf=confs$conf
design=confs$design

system.file("extdata", "Helas3.peak.bed", package="ChIPComp")
setwd("/Library/Frameworks/R.framework/Versions/3.4/Resources/library/ChIPComp/extdata/")
countSet=makeCountSet(conf,design,filetype="bed", species="hg19",binsize=1000)

plot(countSet)

countSet=ChIPComp(countSet)

print(countSet)

# Real data
conf
setwd("/Volumes/FlyPeaks/FlyPeaks")
SampleSheet  <-read.csv("samplesheet/samplesheet_SLX14438_hs_ER_DBA.csv")
SampleSheet<-SampleSheet[c(1,4,3,6,8,9)]
colnames(SampleSheet)<-colnames(conf)
conf<-SampleSheet
design=as.data.frame(lapply(conf[,c("condition","factor")],as.numeric))-1
design=as.data.frame(model.matrix(~condition,design))

#Copied and renamed them sa ChIPComp doesn't like .narrowPeak
#Then used "for bed in *; do cut -f 1-6 $bed > $bed.cut.bed ; done" to convert

conf$peaks<-gsub(conf$peaks,pattern="narrowPeak",replacement="bed.cut.bed")
conf$peaks<-gsub(conf$peaks,pattern="ER/",replacement="ER/bed/")

countSet=makeCountSet(conf,design,filetype="bam", species="hg19")
countSet=ChIPComp(countSet)
countSet$db[c(10:15,17)]

c1<-rowMeans(countSet$db[c(4,6,9)])
c0<-rowMeans(countSet$db[c(5,7,8)])
c1<-c1[]

M<-log(c1*c0)
A<- -log(c1/c0,base=2)

png("plots/049_ChIPComp.png")
par(mar=c(5.1,5.1,4.1,2.1))
plot(M[countSet$db$pvalue.wald>0.05],A[countSet$db$pvalue.wald>0.05],pch=20,
        main="ChIPComp",
        ylab = expression('log'[2]*' Differential ChIP'),
        xlab = expression("log"[10]~"Mean of Normalized Counts"))
points(M[countSet$db$pvalue.wald<0.05],A[countSet$db$pvalue.wald<0.05],pch=20,col="red")
abline(h=0)
dev.off()
