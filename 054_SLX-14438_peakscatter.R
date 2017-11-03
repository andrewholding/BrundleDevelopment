######################
#
# Functions 
#
######################

setwd("/Volumes/FlyPeaks/FlyPeaks")
source('package/brundle.R')

jg.plotScatter<-function(peaks,samples,yaxis,samplesCombinations)
{
    peaks$count <- 1
    peaks_aggregate <- aggregate(count ~ ., peaks, FUN = sum)
    colormap<-rainbow(255,start=0,end=4/6)
    colormap<-rev(colormap)
    peaks_aggregate$Col <- colormap[as.numeric(cut(log(peaks_aggregate$count,2),breaks=255))]
    peaks_aggregate<-peaks_aggregate[order(peaks_aggregate$count),]
    plot(log(peaks_aggregate[,c(1,2)]), col=peaks_aggregate$Col,pch=20, cex=0.3, xlim=c(0,7),ylim=c(0,7),axes=FALSE)
    if (samples > 1) {Axis(side=2, labels=FALSE)} else {
        Axis(side=2, labels=TRUE)
        mtext("log(Counts)", side=2, line=3)
        }
    if (yaxis !=TRUE) {Axis(side=1, labels=FALSE)} else {
        Axis(side=1, labels=TRUE,ylab="Log(Counts)")
        mtext("log(Counts)", side=1, line=3)
    }
    text(1.5,6.5,paste0(samplesCombinations[,samples][1],"-",samplesCombinations[,samples][2]),cex=1.5)
    sampleCor<-cor(log(peaks[1]),log(peaks[2]))
    text(2,6.0,paste0("Correlation = ",signif(sampleCor,3)),cex=1.0)
        abline(a=0,b=1,col="grey")
    return(cor(log(peaks_aggregate[,c(1,2)])))
}


######################
#
# Settings
#
######################

jg.controlMinOverlap      <- 5
jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX14438_hs_CTCF_DBA.csv"
jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX14438_hs_ER_DBA.csv"
jg.treatedCondition       =  "Fulvestrant"
jg.untreatedCondition     =  "none"

######################
#
# Main Code
#
######################

filename<-"Rdata/048_SLX-14438_dba_human_ER_CTCF.rda"
if(!file.exists(filename)){
    dbaExperiment <- jg.getDba(jg.experimentSampleSheet, bRemoveDuplicates=TRUE)
    dbaControl    <- jg.getDba(jg.controlSampleSheet, bRemoveDuplicates=TRUE)
    save(dbaExperiment,dbaControl,file=filename)
} else {
    load(filename)
}

jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
librarySize<-as.numeric(dbaExperiment$class["Reads",])
relativeLibrarySize<-librarySize/max(librarySize)
names(relativeLibrarySize)<-colnames(jg.experimentPeakset[c(-1:-3)])





par(mar=c(0.5, 0.5, 0.2, 0.2), mfrow=c(2,3),
    oma = c(4, 4, 0.2, 0.2))
samplesCombinations<-combn(c("1a","2a","3a"),2)
for (samples in seq(1,length(samplesCombinations)/2))  {
    peaksLibrarySize<-relativeLibrarySize[samplesCombinations[,samples]]
    peaks<-jg.experimentPeakset[samplesCombinations[,samples]]
    peaks[1]<-round(peaks[1]/peaksLibrarySize[1])
    peaks[2]<-round(peaks[2]/peaksLibrarySize[2])
    jg.plotScatter(peaks,samples, F,samplesCombinations)
    }
samplesCombinations<-combn(c("1b","2b","3b"),2)
for (samples in seq(1,length(samplesCombinations)/2))  {
    peaksLibrarySize<-relativeLibrarySize[samplesCombinations[,samples]]
    peaks<-jg.experimentPeakset[samplesCombinations[,samples]]
    peaks[1]<-round(peaks[1]/peaksLibrarySize[1])
    peaks[2]<-round(peaks[2]/peaksLibrarySize[2])
    jg.plotScatter(peaks,samples,T,samplesCombinations)
}

