
######################
#
# Functions 
#
######################

setwd("/Volumes/FlyPeaks/FlyPeaks")
source('package/brundle.R')

######################
#
# Main Code
#
######################

filename<-"Rdata/032_SLX-14229_dba_human_ER_CTCF_normalised.rda"
if(!file.exists(filename)){
  
  dbaSummits                <- 200
  jg.controlMinOverlap      <- 5
  jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX14229_hs_CTCF_DBA.csv"
  jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX14229_hs_ER_DBA.csv"
  jg.treatedCondition       =  "Fulvestrant"
  jg.untreatedCondition     =  "none"
  
  filename_dba<-"Rdata/032_SLX-14229_dba_human_ER_CTCF.rda"
  if(!file.exists(filename)){
    dbaExperiment <- jg.getDba(jg.experimentSampleSheet,dbaSummits)
    dbaControl    <- jg.getDba(jg.controlSampleSheet,   dbaSummits)
    save(dbaExperiment,dbaControl,file=filename)
  } else {
    load(filename)
  }
  
  
  #Condensed code from previous script
  jg.sampleIds <- jg.getSampleIds(jg.controlSampleSheet)
  jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
  jg.controlPeakset    <- jg.dbaGetPeakset(dbaControl)
  jg.controlCountsTreated<-jg.getControlCounts(jg.controlPeakset, 
                                               jg.controlSampleSheet,
                                               jg.treatedCondition)
  jg.controlCountsUntreated<-jg.getControlCounts(jg.controlPeakset,
                                                 jg.controlSampleSheet,
                                                 jg.untreatedCondition)
  jg.untreatedNames <- names(jg.controlCountsUntreated)
  jg.treatedNames   <- names(jg.controlCountsTreated)
  jg.coefficient<-jg.getNormalizationCoefficient(jg.controlCountsTreated,
                                                 jg.controlCountsUntreated)
  jg.correctionFactor<-jg.getCorrectionFactor(jg.experimentSampleSheet,
                                              jg.treatedNames,
                                              jg.untreatedNames)
  jg.experimentPeaksetNormalised<-jd.applyNormalisation(jg.experimentPeakset,
                                                        jg.coefficient,
                                                        jg.correctionFactor,
                                                        jg.treatedNames)
  jg.dba_SLX14229 <- DiffBind:::pv.resetCounts(dbaExperiment,
                                      jg.experimentPeaksetNormalised)

  save(jg.dba_SLX14229,file=filename)
} else {
  load(filename)
}

jg.dba_analysis_SLX14229<-dba.analyze(jg.dba_SLX14229)

png("plots/032_SLX-14229_DiffBind_Analysis.png")
  dba.plotMA(jg.dba_analysis_SLX14229,bFlip=TRUE)
dev.off()


filename<-"Rdata/032_SLX-8047_dba_normalised.rda"
if(!file.exists(filename)){
  
  dbaSummits                <- 200
  jg.controlMinOverlap      <- 5
  jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX8047_dm.csv"
  jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX8047_hs.csv"
  jg.treatedCondition       =  "Fulvestrant"
  jg.untreatedCondition     =  "none"
  
  filename_dba<-"Rdata/032_SLX-8047_dba_hs.rda"
  if(!file.exists(filename)){
    dbaExperiment <- jg.getDba(jg.experimentSampleSheet,dbaSummits)
    dbaControl    <- jg.getDba(jg.controlSampleSheet,   dbaSummits)
    save(dbaExperiment,dbaControl,file=filename)
  } else {
    load(filename)
  }
  
  
  #Condensed code from previous script
  jg.sampleIds <- jg.getSampleIds(jg.controlSampleSheet)
  jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
  jg.controlPeakset    <- jg.dbaGetPeakset(dbaControl)
  jg.controlCountsTreated<-jg.getControlCounts(jg.controlPeakset, 
                                               jg.controlSampleSheet,
                                               jg.treatedCondition)
  jg.controlCountsUntreated<-jg.getControlCounts(jg.controlPeakset,
                                                 jg.controlSampleSheet,
                                                 jg.untreatedCondition)
  jg.untreatedNames <- names(jg.controlCountsUntreated)
  jg.treatedNames   <- names(jg.controlCountsTreated)
  jg.coefficient<-jg.getNormalizationCoefficient(jg.controlCountsTreated,
                                                 jg.controlCountsUntreated)
  jg.correctionFactor<-jg.getCorrectionFactor(jg.experimentSampleSheet,
                                              jg.treatedNames,
                                              jg.untreatedNames)
  jg.experimentPeaksetNormalised<-jd.applyNormalisation(jg.experimentPeakset,
                                                        jg.coefficient,
                                                        jg.correctionFactor,
                                                        jg.treatedNames)
  jg.dba_SLX8047 <- DiffBind:::pv.resetCounts(dbaExperiment,
                                               jg.experimentPeaksetNormalised)
  
  save(jg.dba_SLX8047,file=filename)
} else {
  load(filename)
}


jg.dba_analysis_SLX8047<-dba.analyze(jg.dba_SLX8047)

png("plots/032_SLX-8047_DiffBind_Analysis.png")
dba.plotMA(jg.dba_analysis_SLX8047,bFlip=TRUE)
dev.off()

#Now re-extract the counts and plots xy scatter of similar peaks.

jg.dba_SLX8047
jg.dba_SLX14229

jg.peakset_SLX8047<-dba.peakset(jg.dba_SLX8047, bRetrieve = T, DataType = DBA_DATA_FRAME)
jg.peakset_SLX14229<-dba.peakset(jg.dba_SLX14229, bRetrieve = T, DataType = DBA_DATA_FRAME)


jg.peakset_SLX8047
jg.peakset_SLX14229

names(jg.peakset_SLX14229)[-c(1:3)]<-paste0("SLX_14229_",substr(names(jg.peakset_SLX14229[-c(1:3)]),2,3))
names(jg.peakset_SLX8047)[-c(1:3)]<-paste0("SLX_8047_",substr(names(jg.peakset_SLX8047[-c(1:3)]),2,3))

jg.peakset_combined<-jg.peakset_SLX14229
jg.peakset_combined<-cbind(jg.peakset_combined, matrix(NA,nrow(jg.peakset_combined),ncol(jg.peakset_SLX8047[-c(1:3)])) ) 
colnames(jg.peakset_combined)[-c(1:9)]<-names(jg.peakset_SLX8047[-c(1:3)])

for (peak in rownames(jg.peakset_combined)) {
  matchingPeak<-which.min(abs(
    rowMeans(cbind(jg.peakset_SLX8047$START[jg.peakset_SLX8047["CHR"]==jg.peakset_SLX14229[peak,"CHR"]],
                   jg.peakset_SLX8047$END[jg.peakset_SLX8047["CHR"]==jg.peakset_SLX14229[peak,"CHR"]]))-
      mean(jg.peakset_SLX14229[peak,]$START,jg.peakset_SLX14229[peak,]$END)
  ))
  jg.peakset_combined[peak,][-c(1:9)]<-jg.peakset_SLX8047[matchingPeak,][-c(1:3)]
}

 #Compared log ratios of SLX14229 to SLX8047
jg.peakset_combined
  

  