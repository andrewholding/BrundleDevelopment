
library(devtools)

install_github("andrewholding/Brundle")
library(Brundle)

setwd("/Volumes/FlyPeakCaseStudy/BrundleDevelopment")

#Tidy naming to ER 

jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX15090_CTCF.csv"
jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX15090_H4K12ac.csv"
jg.treatedCondition       =  "Estrogen"
jg.untreatedCondition     =  "none"

#####
#
# Main Code
#
#####

filename<-"Rdata/066_SLX-15090_dba_human_ER_CTCF.rda"
if(!file.exists(filename)){
    dbaExperiment <- jg.getDba(jg.experimentSampleSheet, bRemoveDuplicates = TRUE)
    dbaControl    <- jg.getDba(jg.controlSampleSheet, bRemoveDuplicates = TRUE)
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
png("plots/066__SLX15090_MAplot.png",pointsize=15)
dba.plotMA(jg.dba_analysis,bFlip=TRUE,th=0.05)
dev.off()

report<-dba.report(jg.dba_analysis)
t.test(report$Conc_Estrogen,report$Conc_none,alternative="greater") #p-value = 0.01834

#Check CTCF 
ctcf.dba.analysis<-dba.analyze(dbaControl)
dba.plotMA(ctcf.dba.analysis,bFlip=TRUE) 

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

controlNames<-c("1a","2a","3a") #The ab is flipped in the samplesheet but has
treatmentNames<-c("1b","2b","3b") #no impact on analysis


jg.plotMargins(2,3)

jg.plotPanelsPlot(treatmentNames,
                  controlNames,
                  dbaExperiment
)


png("plots/066_SLX15091_scatter.png", pointsize = 15)
#Worse correlation is 2a-3a (control sample)
controlNames<-c("2a","3a")
jg.plotMargins(1,1)

jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
relativeLibrarySize<-jg.getLibrarySize(dbaExperiment,
                                       jg.experimentPeakset)
peaksLibrarySize<-relativeLibrarySize[controlNames]
peaks<-jg.experimentPeakset[controlNames]
peaks[1]<-round(peaks[1]/peaksLibrarySize[1])
peaks[2]<-round(peaks[2]/peaksLibrarySize[2])
x<-c("Rep2","Rep3")
names(x)<-controlNames
x<-as.matrix(x)
jg.plotScatter(peaks,controlNames, T,t(x),FALSE)
Axis(side=2, labels=TRUE,ylab="Log(Counts)")
mtext("log(Counts)", side=2, line=3)
dev.off()


##Bind at prviously reported sites
dba.report(jg.dba_analysis)
countSites<-dba.report(jg.dba_analysis, th=1)
countSites[countSites$Fold > 0] #3909 sites (Note Fold is flipped this is DOWN)
countSites[countSites$Fold < 0] #9899 sites


report <- dba.report(jg.dba_analysis, th=1,
                     DataType=DBA_DATA_FRAME)
scores <- -10*(log10(report$FDR))
sites  <- cbind(report[,1:3],rownames(report),scores)


#boxplot
par(mfrow=c(1,3))
x<-report[report$Chr=='chr2',] 
x["6687",] #GREB1
barplot(2^as.matrix(x["6687",c("Conc_none","Conc_Estrogen")]), main="GREB1",xlab="Condition",ylab="Read Depth",names.arg=c("Ctrl","E2"))
x<-report[report$Chr=='chr10',] 
x["1295",] #CXCL12
barplot(2^as.matrix(x["1295",c("Conc_none","Conc_Estrogen")]), main="CXCL12",xlab="Condition",ylab="Read Depth",names.arg=c("Ctrl","E2"))
x<-report[report$Chr=='chr22',] 
x["8523",] #XBP1
barplot(2^as.matrix(x["8523",c("Conc_none","Conc_Estrogen")]), main="XBP1",xlab="Condition",ylab="Read Depth",names.arg=c("Ctrl","E2"))



#Not used as quality of plots is only from bed files which doesn't include
#read data.

#source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
#biocLite("ChIPseeker")
#biocLite("clusterProfiler")

#library(ChIPseeker)
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#library(clusterProfiler)
#files<-as.character(read.csv(jg.experimentSampleSheet)$Peaks)

#peak_none <- readPeakFile(files[1])
#tagMatrix_none <- getTagMatrix(peak_none, windows=promoter)
#tagHeatmap(tagMatrix_none, xlim=c(-3000, 3000), color="red")
#plotAvgProf(tagMatrix_none, xlim=c(-3000, 3000), 
#            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")


#peak_E2 <- readPeakFile(files[2])
#tagMatrix_E2 <- getTagMatrix(peak_E2, windows=promoter)
#tagHeatmap(tagMatrix_E2, xlim=c(-3000, 3000), color="red")
#plotAvgProf(tagMatrix_E2, xlim=c(-3000, 3000), 
 #           xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

#tagMatrixList<-list(tagMatrix_none,tagMatrix_E2)
#plotAvgProf( tagMatrixList,xlim=c(-3000, 3000), 
 #            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

#tagHeatmap(tagMatrix_none, xlim=c(-3000, 3000), color="blue")
#tagHeatmap(tagMatrix_E2, xlim=c(-3000, 3000), color="red")


profile<-read.table("SLX-15090/Homer_Profile.txt",sep="\t")
profile_ER<-as.matrix(profile[,c(1,2)])
png("plots/066_HomerTSSProfile.png",pointsize=15)
plot(profile_ER,t="l",col="red", xlab="Distance from TSS",ylab="Read Depth",main="H4K12ac Profile")
profile_Ctrl<-as.matrix(profile[,c(1,5)])
lines(profile_Ctrl,col="blue")
legend("topright",legend = c("Estrogen", "Control"), col=c("red","blue"),pch = c(19,19))
abline(v=0, lty=3)
dev.off()
