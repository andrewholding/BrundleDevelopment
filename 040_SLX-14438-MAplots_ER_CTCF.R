######################
#
# Functions 
#
######################

source('package/brundle.R')

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
png("plots/040_SLX-14438_DiffBind_Analysis.png")
dba.plotMA(jg.dba_analysis,bFlip=TRUE)
dev.off()

#Analyze and plot with Diffbind no correction     
dba_analysis<-dba.analyze(dbaExperiment)
png("plots/040_SLX-14438_DiffBind_Analysis_no_correction.png")
dba.plotMA(dba_analysis,bFlip=TRUE)
dev.off()















#Original vs Normalised data for comparison
png("plots/040_SLX-14438_DiffBind_Analysis_Before_After.png")
par(mfrow=c(1,2))
dba.plotMA(dba_analysis,bFlip=TRUE)
dba.plotMA(jg.dba_analysis,bFlip=TRUE)
par(mfrow=c(1,1))
dev.off()

#FRiP Comparision for manuscript

mean(as.numeric(dbaExperiment$SN))    #2.5%  Reads in ER Peaks 
mean(as.numeric(dbaControl$SN))       #46% Reads in CTCF Peaks

