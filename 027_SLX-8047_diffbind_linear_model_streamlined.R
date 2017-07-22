
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
dbaSummits                <- 200
jg.controlMinOverlap      <- 5
jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX8047_dm.csv"
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

dbaExperiment <- jg.getDba(jg.experimentSampleSheet,dbaSummits)
dbaControl    <- jg.getDba(jg.controlSampleSheet,   dbaSummits)

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
jg.MAplot(jg.experimentPeakset,jg.controlPeakset,jg.untreatedNames,jg.treatedNames,jg.coefficient)

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
dba.plotMA(jg.dba_analysis)

#Original vs Normalised data for comparison
dba_analysis<-dba.analyze(dbaExperiment)
par(mfrow=c(1,2))
dba.plotMA(dba_analysis)
dba.plotMA(jg.dba_analysis)
par(mfrow=c(1,1))

