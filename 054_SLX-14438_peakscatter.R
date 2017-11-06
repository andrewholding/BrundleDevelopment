######################
#
# Functions 
#
######################

setwd("/Volumes/FlyPeaks/FlyPeaks")
source('package/brundle.R')

jg.plotScatter<-function(peaks,samples,yaxis,samplesCombinations,noYLabel=FALSE)
{
    peaks$count <- 1
    peaks_aggregate <- aggregate(count ~ ., peaks, FUN = sum)
    colormap<-rainbow(255,start=0,end=4/6)
    colormap<-rev(colormap)
    peaks_aggregate$Col <- colormap[as.numeric(cut(log(peaks_aggregate$count,2),breaks=255))]
    peaks_aggregate<-peaks_aggregate[order(peaks_aggregate$count),]
    plot(log(peaks_aggregate[,c(1,2)]), col=peaks_aggregate$Col,pch=20, cex=0.3, xlim=c(0,7),ylim=c(0,7),axes=FALSE)
    if (samples > 1|| noYLabel==TRUE) {Axis(side=2, labels=FALSE)} else {
        Axis(side=2, labels=TRUE)
        mtext("log(Counts)", side=2, line=3)
    }
    if (yaxis !=TRUE ) {Axis(side=1, labels=FALSE)} else {
        Axis(side=1, labels=TRUE,ylab="Log(Counts)")
        mtext("log(Counts)", side=1, line=3)
    }
    text(1.5,6.5,paste0(samplesCombinations[,samples][1],"-",samplesCombinations[,samples][2]),cex=1.5)
    sampleCor<-cor(log(peaks[1]),log(peaks[2]))
    text(2,6.0,paste0("Correlation = ",signif(sampleCor,3)),cex=1.0)
        abline(a=0,b=1,col="grey")
    return(cor(log(peaks_aggregate[,c(1,2)])))
}


jg.plotMargins<-function(a,b)
{
    par(mar=c(0.5, 0.5, 0.2, 0.2), mfrow=c(a,b),
        oma = c(4, 4, 0.2, 0.2))
}
jg.plotPanelsPlot<-function(treatmentNames,controlNames,dbaExperiment,noYLabel=FALSE)
{
    
    jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
    relativeLibrarySize<-jg.getLibrarySize(dbaExperiment,
                                            jg.experimentPeakset)
    
    samplesCombinations<-combn(controlNames,2)
    for (samples in seq(1,length(samplesCombinations)/2))  {
        peaksLibrarySize<-relativeLibrarySize[samplesCombinations[,samples]]
        peaks<-jg.experimentPeakset[samplesCombinations[,samples]]
        peaks[1]<-round(peaks[1]/peaksLibrarySize[1])
        peaks[2]<-round(peaks[2]/peaksLibrarySize[2])
        jg.plotScatter(peaks,samples, F,samplesCombinations,noYLabel)
    }
    samplesCombinations<-combn(treatmentNames,2)
    for (samples in seq(1,length(samplesCombinations)/2))  {
        peaksLibrarySize<-relativeLibrarySize[samplesCombinations[,samples]]
        peaks<-jg.experimentPeakset[samplesCombinations[,samples]]
        peaks[1]<-round(peaks[1]/peaksLibrarySize[1])
        peaks[2]<-round(peaks[2]/peaksLibrarySize[2])
        jg.plotScatter(peaks,samples,T,samplesCombinations,noYLabel)
    }
 }


jg.getLibrarySize<-function(dbaExperiment,jg.experimentPeakset) {
    librarySize<-as.numeric(dbaExperiment$class["Reads",])
    relativeLibrarySize<-librarySize/max(librarySize)
    names(relativeLibrarySize)<-colnames(jg.experimentPeakset[c(-1:-3)])
    return(relativeLibrarySize)
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


### SLX14438

filename<-"Rdata/048_SLX-14438_dba_human_ER_CTCF.rda"
load(filename)

controlNames<-c("1a","2a","3a")
treatmentNames<-c("1b","2b","3b")

jg.plotMargins(2,3)

jg.plotPanelsPlot(treatmentNames,
                  controlNames,
                  dbaExperiment
                  )


### SLX8047
filename<-"Rdata/026_SLX-8047_dba.counts_human.rda"
load(filename)
dbaExperiment_SLX8047<-dba

controlNames_SLX8047<-c("1b","2b","3b","4b")
treatmentNames_SLX8047<-c("1a","2a","3a","4a")

jg.plotMargins(2,6)

jg.plotPanelsPlot(treatmentNames_SLX8047,
                  controlNames_SLX8047,
                  dbaExperiment_SLX8047
)

### SLX-12998

filename<-"Rdata/047_SLX-12998_dba_human.rda"
load(filename)
dbaExperiment_SLX12998<-dbaFullExperiment

jg.plotMargins(2,6)

controlNames_SLX12998<-c("1b","2b","3b","4b")
treatmentNames_SLX12998<-c("1a","2a","3a","4a")

jg.plotPanelsPlot(treatmentNames_SLX12998,
                  controlNames_SLX12998,
                  dbaExperiment_SLX12998
                )
#### Figure for paper

jg.plotMargins(2,3)
par(mfcol=c(2,3))

#SLX-14438
#2a-3a and 1b-3b
controlNames<-c("2a","3a")
treatmentNames<-c("1b","3b")
jg.plotPanelsPlot(treatmentNames,
                  controlNames,
                  dbaExperiment)

#SLX-8047
#3b-4b and 1a-4a
controlNames<-c("3b","4b")
treatmentNames<-c("1a","4a")
jg.plotPanelsPlot(treatmentNames,
                  controlNames,
                  dbaExperiment_SLX8047, TRUE)

#SLX-12998
#3b-4b and 1a-4a
controlNames<-c("3b","4b")
treatmentNames<-c("1a","4a")
jg.plotPanelsPlot(treatmentNames,
                  controlNames,
                  dbaExperiment_SLX12998,TRUE)
