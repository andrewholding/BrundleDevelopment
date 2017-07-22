
######################
#
# Functions 
#
######################

source("package/brundle.R")

######################
#
# Settings
#
######################

setwd("/Volumes/FlyPeaks/FlyPeaks")
dbaSummits                <- 200
jg.controlMinOverlap      <- 5
jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX8047_hs.csv"
jg.treatedCondition       =  "Fulvestrant"
jg.untreatedCondition     =  "none"

######################
#
# Main Code
#
######################

# To prepare data for this code, first align the sequencing data to a combined
# genome of the experimental and control species chromatin. Then split the bam
# agliments by species, reindex and call peaks.
#
# The workflow that follows takes the experimental BAM files and peaks from a
# standard DiffBind samplesheet along with a matching samplessheet with the
# approriate BAM and bedfiles for the control peaks.
#
# From these the script generates a normalisation factor and then returns the
# to output as a the DiffBind object for downstream processing.
#

#Load experimental data
dbaExperiment <- jg.getDba(jg.experimentSampleSheet,dbaSummits)

#Load Sample Ids from control sample sheet.
jg.sampleIds <- jg.getSampleIds(jg.experimentSampleSheet)

#We want a list of samplesheets to try for our control regions
jg.controlSampleSheet<-c()
for (jg.titration in c("A1","A2","A5","A10","A15", "A20","A30","A40","A50")) {
  
  jg.controlSampleSheet[jg.titration]     <- paste0("samplesheet/titration/samplesheet_SLX8047_dm_",jg.titration,".csv")
     
  dbaControl    <- jg.getDba(jg.controlSampleSheet[jg.titration],   dbaSummits)
  
  ## Extract Peak set from DiffBind
  jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
  jg.controlPeakset    <- jg.dbaGetPeakset(dbaControl)
  
  
  #Get counts for each condition
  jg.controlCountsTreated<-jg.getControlCounts(jg.controlPeakset, 
                                               jg.controlSampleSheet[jg.titration],
                                               jg.treatedCondition)
  jg.controlCountsUntreated<-jg.getControlCounts(jg.controlPeakset,
                                                 jg.controlSampleSheet[jg.titration],
                                                 jg.untreatedCondition)
  
  #Get sample names for conditions
  jg.untreatedNames <- names(jg.controlCountsUntreated)
  jg.treatedNames   <- names(jg.controlCountsTreated)
  
  ##Plot showing normalization calculation (Optional)
  png (paste0("plots/031_DBA_normalization_plot_",jg.titration,".png"))
    jg.plotNormalization(jg.controlCountsTreated,
                         jg.controlCountsUntreated)
  dev.off()
  ##Get Normalization Coefficient
  jg.coefficient<-jg.getNormalizationCoefficient(jg.controlCountsTreated,
                                                 jg.controlCountsUntreated)
  
  
  #Check by MA plot (Optional)
  png (paste0("plots/031_MA_corrected_",jg.titration,".png"))
   jg.plotMA(jg.experimentPeakset,jg.controlPeakset,jg.untreatedNames,jg.treatedNames,jg.coefficient)
  dev.off()
  #DeSeq called by Diffbind will divide though by library size to normalise, 
  #and therefore partially undo our work. Solution is to correct our normalisation
  #factor to with this correction factor to remove library size part of our 
  #NormalisationFactor.
  
  jg.correctionFactor<-jg.getCorrectionFactor(jg.experimentSampleSheet,
                                              jg.treatedNames,
                                              jg.untreatedNames
                                              )
  
  #Apply coefficent and control factor
  jg.experimentPeaksetNormalised<-jd.applyNormalisation(jg.experimentPeakset,
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
  png (paste0("plots/031_DBA_corrected_",jg.titration,".png"))
    dba.plotMA(jg.dba_analysis)
  dev.off()
} 
#Unnormalised data for comparison
dba_analysis<-dba.analyze(dbaExperiment)
png (paste0("plots/031_DBA_uncorrected.png"))
  dba.plotMA(dba_analysis)
dev.off()
