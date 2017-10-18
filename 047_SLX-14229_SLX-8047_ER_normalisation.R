
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

dbaSummits                <- 200
jg.controlMinOverlap      <- 5
jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX14438_hs_CTCF_safe.csv"
jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX14438_hs_ER_DBA.csv"
jg.treatedCondition       =  "Fulvestrant"
jg.untreatedCondition     =  "none"

######################
#
# Main Code
#
######################

filename<-"Rdata/047_SLX-14438_dba_human_ER_CTCF.rda"
if(!file.exists(filename)){
  dbaExperiment <- jg.getDba(jg.experimentSampleSheet,dbaSummits)
  dbaControl    <- jg.getDba(jg.controlSampleSheet,   dbaSummits)
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

##Plot showing normalization calculation (Optional)
  jg.plotNormalization(jg.controlCountsTreated,
                     jg.controlCountsUntreated)

 ##Get Normalization Coefficient
jg.coefficient<-jg.getNormalizationCoefficient(jg.controlCountsTreated,
                                               jg.controlCountsUntreated)


#Check by MA plot (Optional)
jg.plotMA(jg.controlPeakset,jg.experimentPeakset,jg.untreatedNames,jg.treatedNames,jg.coefficient)


#As we're going to normalise a second data set we don't need to use a DiffBind correction factor
jg.correctionFactor<-1


#Apply coefficent and control factor
jg.experimentPeaksetNormalised<-jg.applyNormalisation(jg.experimentPeakset,
                                                      jg.coefficient,
                                                      jg.correctionFactor,
                                                      jg.treatedNames
)


#Return values to Diffbind object
jg.dba <- DiffBind:::pv.resetCounts(dbaExperiment,
                                    jg.experimentPeaksetNormalised
)

#jg.dba is now normalised a list of ER peaks and their expected ratio


#####
#
#Rerun the normalisation code to get the ER ratio
#
######


jg.normalisedExperimentPeakset <- jg.dbaGetPeakset(jg.dba)


#Get counts for each condition
jg.ERCountsTreated<-jg.getControlCounts(jg.normalisedExperimentPeakset, 
                                             jg.experimentSampleSheet,
                                             jg.treatedCondition)
jg.ERCountsUntreated<-jg.getControlCounts(jg.normalisedExperimentPeakset,
                                          jg.experimentSampleSheet,
                                               jg.untreatedCondition)

#Get sample names for conditions
jg.untreatedNames <- names(jg.ERCountsTreated)
jg.treatedNames   <- names(jg.ERCountsUntreated)

jg.plotMA(jg.normalisedExperimentPeakset,jg.normalisedExperimentPeakset,jg.untreatedNames,jg.treatedNames,1)


##Establish normalised ER ratio as we know it isn't 1

jg.plotNormalization(jg.ERCountsTreated,
                     jg.ERCountsUntreated)

jg.ERratio<-jg.getNormalizationCoefficient(jg.ERCountsTreated,
                                               jg.ERCountsUntreated)

########
#
# Apply to ER only data
#
########

jg.EROnlyExperimentSampleSheet <- "samplesheet/samplesheet_SLX8047_hs.csv"

filename<-"Rdata/047_SLX-SLX-8047_dba_human_ER.rda"
if(!file.exists(filename)){
    dbaERonlyExperiment <- jg.getDba(jg.EROnlyExperimentSampleSheet,dbaSummits)
    save(dbaERonlyExperiment,file=filename)
} else {
    load(filename)
}

#Get peakset

jg.EROnlyExperimentPeakset <- jg.dbaGetPeakset(dbaERonlyExperiment)


#Get counts for each condition
jg.EROnlyCountsTreated<-jg.getControlCounts(jg.EROnlyExperimentPeakset, 
                                            jg.EROnlyExperimentSampleSheet,
                                        jg.treatedCondition)
jg.EROnlyCountsUntreated<-jg.getControlCounts(jg.EROnlyExperimentPeakset,
                                              jg.EROnlyExperimentSampleSheet,
                                          jg.untreatedCondition)

jg.plotNormalization(jg.EROnlyCountsTreated,
                     jg.EROnlyCountsUntreated)

#Get sample names for conditions
jg.untreatedNames <- names(jg.EROnlyCountsTreated)
jg.treatedNames   <- names(jg.EROnlyCountsUntreated)



jg.EROnlyExperimentRatio<-jg.getNormalizationCoefficient(jg.EROnlyCountsTreated,
                                           jg.EROnlyCountsUntreated)

#Create new ratio 
jg.ERNormalisationRatio<-jg.ERratio/jg.EROnlyExperimentRatio

#PlotMA (Optional)

jg.plotMA(jg.EROnlyExperimentPeakset,jg.EROnlyExperimentPeakset,jg.treatedNames,jg.untreatedNames,1)
jg.plotMA(jg.EROnlyExperimentPeakset,jg.EROnlyExperimentPeakset,jg.treatedNames,jg.untreatedNames,1/jg.ERNormalisationRatio)

