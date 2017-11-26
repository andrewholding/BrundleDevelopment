
library(devtools)

install_github("andrewholding/Brundle")
library(Brundle)

setwd("/Volumes/FlyPeakCaseStudy/BrundleDevelopment")

jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX15091_CTCF.csv"
jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX15091_ER.csv"
jg.treatedCondition       =  "Estrogen"
jg.untreatedCondition     =  "none"

#####
#
# Main Code
#
#####

filename<-"Rdata/065_SLX-15091_dba_human_ER_CTCF.rda"
if(!file.exists(filename)){
    dbaExperiment <- jg.getDba(jg.experimentSampleSheet)
    dbaControl    <- jg.getDba(jg.controlSampleSheet)
    save(dbaExperiment,dbaControl,file=filename)
} else {
    load(filename)
}

jg.dba<-Brundle(
            dbaExperiment,
            dbaControl,
            jg.treatedCondition,
            jg.untreatedCondition,
            jg.experimentSampleSheet,
            jg.controlSampleSheet )
    
jg.dba_analysis<-dba.analyze(jg.dba)
png("plots/066_SLX15091_MAplot.png",pointsize=15)
dba.plotMA(jg.dba_analysis,bFlip=TRUE,th=0.05)
dev.off()

#Check CTCF 
ctcf.dba.analysis<-dba.analyze(dbaControl)
dba.plotMA(ctcf.dba.analysis,bFlip=TRUE) 

dba.report(jg.dba_analysis)


#IDR for manuscript


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

controlNames<-c("1a","2a","3b") #The ab is flipped in the samplesheet but has
treatmentNames<-c("1b","2b","3a") #no impact on analysis


jg.plotMargins(2,3)

jg.plotPanelsPlot(treatmentNames,
                  controlNames,
                  dbaExperiment
)


png("plots/066_SLX15091_scatter.png", pointsize = 15)
#Worse correlation is 1a-3b (control sample)
controlNames<-c("1a","3b")
jg.plotMargins(1,1)

jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
relativeLibrarySize<-jg.getLibrarySize(dbaExperiment,
                                       jg.experimentPeakset)
peaksLibrarySize<-relativeLibrarySize[controlNames]
peaks<-jg.experimentPeakset[controlNames]
peaks[1]<-round(peaks[1]/peaksLibrarySize[1])
peaks[2]<-round(peaks[2]/peaksLibrarySize[2])
x<-c("Rep1","Rep3")
names(x)<-controlNames
x<-as.matrix(x)
jg.plotScatter(peaks,controlNames, T,t(x),FALSE)
Axis(side=2, labels=TRUE,ylab="Log(Counts)")
mtext("log(Counts)", side=2, line=3)
dev.off()


###Export for Homer
dba.report(jg.dba_analysis)


report <- dba.report(jg.dba_analysis, th=.05,
                     DataType=DBA_DATA_FRAME)
scores <- -10*(log10(report$FDR))
sites  <- cbind(report[,1:3],rownames(report),scores)

#Remember the data is flipped by DiffBind
gains  <- report$Fold < 0
losses <- report$Fold > 0

write.table(sites[gains,],
               "bed/ER_45minsGains.bed", quote=FALSE, sep="\t",
               row.names=FALSE, col.names=FALSE)

dev.off()
#boxplot
par(mfrow=c(1,3))
x<-report[report$Chr=='chr17',] 
x["9769",] #RARa
barplot(2^as.matrix(x["9769",c("Conc_none","Conc_Estrogen")]), main="RARa",xlab="Condition",ylab="Read Depth",names.arg=c("Ctrl","E2"))
x<-report[report$Chr=='chr21',] 
x["14816",] #NRIP
barplot(2^as.matrix(x["14816",c("Conc_none","Conc_Estrogen")]), main="NRIP",xlab="Condition",ylab="Read Depth",names.arg=c("Ctrl","E2"))
x<-report[report$Chr=='chr22',] 
x["15381",] #XBP1
barplot(2^as.matrix(x["15381",c("Conc_none","Conc_Estrogen")]), main="XBP1",xlab="Condition",ylab="Read Depth",names.arg=c("Ctrl","E2"))




