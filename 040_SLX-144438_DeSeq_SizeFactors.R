
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
png("plots/040_SLX-14438_DeSeq_Analysis_ER.png")
  jg.plotDeSeq(jg.controlResultsDeseq,
               p=0.01, 
               title.main="Fold-change in CTCF binding",
               flip=T
               )
dev.off()
png("plots/040_SLX-14438_DeSeq_Analysis_CTCF.png")
jg.plotDeSeq(jg.experimentResultsDeseq,
             p=0.01,
             title.main="Fold-change in ER binding",
             flip=T
             )
dev.off()

#Repeat not using out control
jg.experimentDeSeqInternal<-jg.runDeSeq(jg.experimentPeaksetDeSeq, jg.conditions, NULL)
jg.experimentResultsDeseqInternal   = results(jg.experimentDeSeqInternal)

png("plots/040_SLX-14438_DeSeq_Analysis_ER_nocorrection.png")
  jg.plotDeSeq(jg.experimentResultsDeseqInternal,
               title.main="Fold-change in ER binding (no correction)",
               p=0.01,
               flip=T
               )
dev.off()


#Draw Combined figure.
png("plots/040_SLX-14438_DeSeq_Analysis_ER_CTCF.png")
  jg.plotDeSeqCombined(jg.controlResultsDeseq,
                       jg.experimentResultsDeseq,
                       title.main="ER and CTCF Binding Folding changes on ER treatment",
                       p=0.01,flip=TRUE)
dev.off()


#Correct in DiffBind
dbaExperimentCorrected<-jg.correctDBASizeFactors(dbaExperiment,
                                                 jg.controlSizeFactors
                                                 )

#Plot DiffBind MA - Before and After
png("plots/040_SLX-14438_DiffBind_DeSeq_Analysis_ER_nocorrection.png")
  par(mar=c(5.1,5.1,4.1,4.1))
  dbaExperimentAnalysis<-dba.analyze(dbaExperiment)
  dba.plotMA(dbaExperimentAnalysis,bFlip=TRUE,cex.lab=1.5, cex.axis=1.5, cex.main=1.25, cex.sub=1.5)
dev.off()
  
png("plots/040_SLX-14438_DiffBind_DeSeq_Analysis_ER.png")
  par(mar=c(5.1,5.1,4.1,4.1))
  dbaExperimentAnalysisCorrected<-dba.analyze(dbaExperimentCorrected)
  dba.plotMA(dbaExperimentAnalysisCorrected,bFlip=TRUE,cex.lab=1.5, cex.axis=1.5, cex.main=1.25, cex.sub=1.5)
dev.off()


######################
#
# Code snippets used to calculate figures for manuscript
#
######################


#Count number of differntially bound CTCF sites to see if it is less than FDR.
nrow(jg.controlResultsDeseq)

#Total number of signficatly changed sites
jg.DeSeqSignificant <- jg.controlResultsDeseq$"padj" < 0.01 & !is.na(jg.controlResultsDeseq$"padj")
nrow(jg.controlResultsDeseq[jg.DeSeqSignificant,])

#Up Only (not < 0 as we have to flip data to plot)
jg.DeSeqSignificantUp <- jg.controlResultsDeseq[jg.DeSeqSignificant,]$"log2FoldChange" < 0
nrow(jg.controlResultsDeseq[jg.DeSeqSignificant,][jg.DeSeqSignificantUp,])

#Down Only
jg.DeSeqSignificantDown <- jg.controlResultsDeseq[jg.DeSeqSignificant,]$"log2FoldChange" > 0
nrow(jg.controlResultsDeseq[jg.DeSeqSignificant,][jg.DeSeqSignificantDown,])


