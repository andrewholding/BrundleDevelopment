setwd("/Volumes/FlyPeaks/FlyPeaks")
source('package/brundle.R')
jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX14438_hs_CTCF_DBA_noER.csv"
jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX14438_hs_ER_DBA.csv"
jg.treatedCondition       =  "Fulvestrant"
jg.untreatedCondition     =  "none"
dbaSummits                <- 200

filename<-"Rdata/059_SLX-14438_dba_human_ER_CTCF_noER.rda"
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

##Get Normalization Coefficient
jg.coefficient<-jg.getNormalizationCoefficient(jg.controlCountsTreated,
                                               jg.controlCountsUntreated)


#Check by MA plot (Optional)
jg.plotMA(jg.experimentPeakset,jg.controlPeakset,jg.untreatedNames,jg.treatedNames,jg.coefficient)

#DeSeq called by Diffbind will divide though by library size to normalise, 
#and therefore partially undo our work. Solution is to correct our normalisation
#factor to with this correction factor to remove library size part of our 
#NormalisationFactor.

jg.correctionFactor<-jg.getCorrectionFactor(jg.experimentSampleSheet,
                                            jg.treatedNames,
                                            jg.untreatedNames
)


#Apply coefficent and control factor
jg.experimentPeaksetNormalised<-jg.applyNormalisation(jg.experimentPeakset,
                                                      jg.coefficient,
                                                      jg.correctionFactor,
                                                      jg.treatedNames
)


#Return values to Diffbind and plot normalised result.
jg.dba <- DiffBind:::pv.resetCounts(dbaExperiment,
                                    jg.experimentPeaksetNormalised
)


#Analyze and plot with Diffbind      
jg.dba_analysis<-dba.analyze(jg.dba)

dba.plotMA(jg.dba_analysis,bFlip=TRUE)
           
#Check CTCF 
ctcf.dba.analysis<-dba.analyze(dbaControl)
 dba.plotMA(ctcf.dba.analysis,bFlip=TRUE)          
 