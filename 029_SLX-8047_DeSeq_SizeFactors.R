
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
# From these the script generates DeSeq sizes factors form the control samples
# and suppiles them to DeSeq as for downstream processing.
#

#The key function is jg.getDba, the saving of file is to speed up the 
#loading of data which can be slow.

filename<-"Rdata/029_SLX-8047_dba_human_drosophila.rda"
if(!file.exists(filename)){
  dbaExperiment <- jg.getDba(jg.experimentSampleSheet,   dbaSummits)
  dbaControl    <- jg.getDba(jg.controlSampleSheet,   dbaSummits)
  save(dbaExperiment,dbaControl,file=filename)
} else {
  load(filename)
}


## Extract Peak set from DiffBind
jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
jg.controlPeakset    <- jg.dbaGetPeakset(dbaControl)

#Convert Peakset to DeSeq Workflow
jg.controlPeaksetDeSeq<-jg.convertPeakset(jg.controlPeakset)

#Establish size factors directly from Control data
jg.controlSizeFactors = estimateSizeFactorsForMatrix(jg.controlPeaksetDeSeq)

#Get conditions dataframe for DeSeq
jg.conditions <- read.csv(file=jg.controlSampleSheet, header=TRUE, sep=",")['Condition']

#Run DeSeq on control
jg.controlDeSeq<-jg.runDeSeq(jg.controlPeaksetDeSeq, jg.conditions,jg.controlSizeFactors)
jg.controlResultsDeseq   = results(jg.controlDeSeq)

#Repeat for experimental conditions

#Convert experiment Peakset to DeSeq Workflow
jg.experimentPeaksetDeSeq<-jg.convertPeakset(jg.experimentPeakset)

#Run DeSeq on experiment
jg.experimentDeSeq<-jg.runDeSeq(jg.experimentPeaksetDeSeq, jg.conditions,jg.controlSizeFactors)
jg.experimentResultsDeseq   = results(jg.experimentDeSeq)

#Plot results
png("plots/029_SLX-8047_DeSeq_Analysis_H2av.png")
  jg.plotDeSeq(jg.controlResultsDeseq,
               p=0.01, 
               title.main="Fold-change in H2av binding"
               )
dev.off()
png("plots/029_SLX-8047_DeSeq_Analysis_ER.png")
jg.plotDeSeq(jg.experimentResultsDeseq,
             p=0.01,
             title.main="Fold-change in ER binding"
             )
dev.off()

#Repeat not using out control
jg.experimentDeSeqInternal<-jg.runDeSeq(jg.experimentPeaksetDeSeq, jg.conditions, NULL)
jg.experimentResultsDeseqInternal   = results(jg.experimentDeSeqInternal)

png("plots/029_SLX-8047_DeSeq_Analysis_ER_nocorrection.png")
  jg.plotDeSeq(jg.experimentResultsDeseqInternal,
               title.main="Fold-change in ER binding (no correction)",
               p=0.01)
dev.off()


#Draw Combined figure.
png("plots/029_SLX-8047_DeSeq_Analysis_HsDm.png")
  jg.plotDeSeqCombined(jg.controlResultsDeseq,
                       jg.experimentResultsDeseq,
                       title.main="ER and H2av Binding Folding changes on ER treatment",
                       p=0.01, flip=T)
dev.off()


#Correct in DiffBind
dbaExperimentCorrected<-jg.correctDBASizeFactors(dbaExperiment,
                                                 jg.controlSizeFactors
                                                 )

#Plot DiffBind MA - Before and After
png("plots/029_SLX-8047_DiffBind_Analysis_ER_nocorrection.png")
  dbaExperimentAnalysis<-dba.analyze(dbaExperiment)
  dba.plotMA(dbaExperimentAnalysis)
dev.off()
  
png("plots/029_SLX-8047_DiffBind_Analysis_ER.png")
  dbaExperimentAnalysisCorrected<-dba.analyze(dbaExperimentCorrected)
  dba.plotMA(dbaExperimentAnalysisCorrected)
dev.off()
