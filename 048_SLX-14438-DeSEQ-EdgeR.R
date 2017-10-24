######################
#
# Functions 
#
######################

setwd("/Volumes/FlyPeaks/FlyPeaks")
source('package/brundle.R')


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




#Load Sample Ids from control sample sheet.
jg.sampleIds <- jg.getSampleIds(jg.controlSampleSheet)

## Extract Peak set from DiffBind
jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
jg.controlPeakset    <- jg.dbaGetPeakset(dbaControl)


#Get counts for each condition
jg.controlCountsTreated<-jg.getControlCounts(jg.controlPeakset, 
                                             jg.controlSampleSheet,
                                             jg.treatedCondition)
jg.controlCountsUntreated<-jg.getControlCounts(jg.controlPeakset,
                                               jg.controlSampleSheet,
                                               jg.untreatedCondition)

#Get sample names for conditions
jg.untreatedNames <- names(jg.controlCountsUntreated)
jg.treatedNames   <- names(jg.controlCountsTreated)


jg.coefficient<-jg.getNormalizationCoefficient(jg.controlCountsTreated,
                                               jg.controlCountsUntreated)
jg.correctionFactor<-jg.getCorrectionFactor(jg.experimentSampleSheet,
                                            jg.treatedNames,
                                            jg.untreatedNames)

#Apply coefficent and control factor
jg.experimentPeaksetNormalised<-jg.applyNormalisation(jg.experimentPeakset,
                                                      jg.coefficient,
                                                      jg.correctionFactor,
                                                      jg.treatedNames)
#Return values to Diffbind and plot normalised result.
jg.dba <- DiffBind:::pv.resetCounts(dbaExperiment,
                                    jg.experimentPeaksetNormalised)


#Analyze and plot with Diffbind      
jg.dba_analysis<-dba.analyze(jg.dba)
par(mar=c(5.1,5.1,4.1,4.1))
dba.plotMA(jg.dba_analysis,bFlip=TRUE, cex.lab=1.5, cex.axis=1.5, cex.main=1.25, cex.sub=1.5)

#Analyze and plot with Diffbind no correction     
dba_analysis_DeSEQ2<-dba.analyze(dbaExperiment,method=DBA_DESEQ2, ,bFullLibrarySize=FALSE)
dba_analysis_EDGER<-dba.analyze(dbaExperiment,method=DBA_EDGER)
dba_analysis_DiffBind<-dba.analyze(dbaExperiment,method=DBA_DESEQ2, ,bFullLibrarySize=TRUE)
dba.plotMA(dba_analysis,bFlip=TRUE,cex.lab=1.5, cex.axis=1.5, cex.main=1.25, cex.sub=1.5)

#Original vs Normalised data for comparison
par(mfrow=c(1,4))
dba.plotMA(dba_analysis_EDGER,bFlip=TRUE, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, method=DBA_EDGER)
dba.plotMA(dba_analysis_DeSEQ2,bFlip=TRUE, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, method=DBA_DESEQ2)
dba.plotMA(dba_analysis_DiffBind,bFlip=TRUE, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, method=DBA_DESEQ2)
dba.plotMA(jg.dba_analysis,bFlip=TRUE,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
par(mfrow=c(1,1))

#Comparison of Reads in Peaks between data sets
boxplot(list(rowMeans(jg.experimentPeakset[-c(1:3)]),rowMeans(jg.controlPeakset[-c(1:3)])))
boxplot(list(rowMeans(jg.experimentPeakset[c(4,6,8)]),rowMeans(jg.controlPeakset[-c(1:3)])))
boxplot(list(rowMeans(jg.experimentPeakset[c(5,6,7)]),rowMeans(jg.controlPeakset[-c(1:3)])))



plotForComparision<-function(dba, MethodName, dbatype,ylimits) {
    df<-as.data.frame(dba.report(dba, method=dbatype, th=1))    
    countSet<-df[as.character(sort(as.numeric(rownames(df)))),]
    peakSet<-jg.dbaGetPeakset(dba)   
    c1<-rowMeans(peakSet[c(4,6,9)])
    c0<-rowMeans(peakSet[c(5,7,8)])
    c1<-c1[]
        

M<-log(c1*c0)
A<- -log(c1/c0,base=2)

A<- -countSet$Fold
#c0<-countSet$Conc_Fulvestrant
#c1<-c1[]


png(paste0("plots/048_",MethodName,".png"), pointsize=20)
par(mar=c(5.1,5.1,4.1,2.1))
plot(M[countSet$FDR>0.05],A[countSet$FDR>0.05],pch=20,
     main=MethodName,
     ylab = expression('log'[2]*' Differential ChIP'),
     xlab = expression("log"[10]~"Mean of Normalized Counts"),cex=0.05,
     xlim=c(0,15),ylim=ylimits
     )
points(M[countSet$FDR<0.05],A[countSet$FDR<0.05],pch=20,col="red",cex=0.05)
abline(h=0)
dev.off()
}

plotForComparision(dba_analysis_EDGER, "EdgeR", DBA_EDGER,c(-3,3))
plotForComparision(dba_analysis_DeSEQ2, "DeSEQ2",DBA_DESEQ2,c(-3,3))
plotForComparision(dba_analysis_DiffBind, "DiffBind",DBA_DESEQ2,c(-6,2))
plotForComparision(jg.dba_analysis, "Brundle",DBA_DESEQ2,c(-6,2))

