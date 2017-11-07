######################
#
# Functions 
#
######################

source('package/brundle.R')

## Modifying plot for publication

jg.plotNormalization <- function(jg.controlCountsTreated,jg.controlCountsUntreated) {
    plot(rowMeans(jg.controlCountsTreated),rowMeans(jg.controlCountsUntreated),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, pch=20,
         xlab="Counts in peak after treatment" ,  ylab="Counts in peak before treatment" ,
         main="Comparision of Counts in peaks",)
    lm1<-lm(rowMeans(jg.controlCountsUntreated) ~ 0 + rowMeans(jg.controlCountsTreated))
    
    abline(c(0,lm1$coef),col="red3")
    print(lm1$coefficients)
    angularcoeff<-lm1$coef[1]
    
    points(rowMeans(jg.controlCountsTreated)*angularcoeff,rowMeans(jg.controlCountsUntreated),pch=20, col="royalblue3" )
    treatment_fit<-rowMeans(jg.controlCountsTreated)*angularcoeff
    lm1<-lm(treatment_fit ~ 0 + rowMeans(jg.controlCountsUntreated))
    abline(c(0,lm1$coef),col="purple")
    legend("topleft",legend=c("Raw", "Normalised"),pch=20,col=c("black","royalblue3"))
}

######################
#
# Settings
#
######################

setwd("/Volumes/FlyPeaks/FlyPeaks")
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

filename<-"Rdata/040_SLX-14438_dba_human_ER_CTCF.rda"
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

##Get Normalization Coefficient
png("plots/040_SLX-14438_Normalization_c.png")
par(mar=c(5.1,5.1,4.1,2.1))
    jg.plotNormalization(jg.controlCountsTreated, jg.controlCountsUntreated)
dev.off()

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
dba.SLX14438_report<-dba.report(jg.dba_analysis, th=1)
filename<-"Rdata/040_SLX-14438_dba_human_ER_CTCF_normalised_report.rda"
save(dba.SLX14438_report,file=filename)


png("plots/040_SLX-14438_DiffBind_Analysis.png")
par(mar=c(5.1,5.1,4.1,4.1))
dba.plotMA(jg.dba_analysis,bFlip=TRUE, cex.lab=1.5, cex.axis=1.5, cex.main=1.25, cex.sub=1.5)
dev.off()

#Analyze and plot with Diffbind no correction     
dba_analysis<-dba.analyze(dbaExperiment)
png("plots/040_SLX-14438_DiffBind_Analysis_no_correction.png")
par(mar=c(5.1,5.1,4.1,4.1))
dba.plotMA(dba_analysis,bFlip=TRUE,cex.lab=1.5, cex.axis=1.5, cex.main=1.25, cex.sub=1.5)
dev.off()

#Original vs Normalised data for comparison
png("plots/040_SLX-14438_DiffBind_Analysis_Before_After.png")
par(mfrow=c(1,2))
dba.plotMA(dba_analysis,bFlip=TRUE, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
dba.plotMA(jg.dba_analysis,bFlip=TRUE,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
par(mfrow=c(1,1))
dev.off()

#Comparison of Reads in Peaks between data sets
boxplot(list(rowMeans(jg.experimentPeakset[-c(1:3)]),rowMeans(jg.controlPeakset[-c(1:3)])))
boxplot(list(rowMeans(jg.experimentPeakset[c(4,6,8)]),rowMeans(jg.controlPeakset[-c(1:3)])))
boxplot(list(rowMeans(jg.experimentPeakset[c(5,6,7)]),rowMeans(jg.controlPeakset[-c(1:3)])))
