
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
jg.experimentResultsDeseq_SLX8047   = results(jg.experimentDeSeq)


png("plots/029_SLX-8047_SLX-14229_DeSeq_Analysis_ER.png")
jg.plotDeSeq(jg.experimentResultsDeseq_SLX8047,
             p=0.01,
             title.main="Fold-change in ER binding",
             flip=T
             )
dev.off()

#Now run for SLX-14229

jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX14229_hs_CTCF_DBA.csv"
jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX14229_hs_ER_DBA.csv"


filename<-"Rdata/028_SLX-14229_dba_human_ER_CTCF.rda"
#As there seems to be seed involved in counting which gives
#a small variance to the the results which is best avoided.

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
jg.experimentResultsDeseq_SLX14229   = results(jg.experimentDeSeq)

#Plot results

png("plots/029_SLX-8047-SLX-14229_DeSeq_Analysis_CTCF.png")
jg.plotDeSeq(jg.experimentResultsDeseq_SLX14229,
             p=0.01,
             title.main="Fold-change in ER binding",
             flip=T
)
dev.off()

library(latticeExtra)

ma.df<-jg.experimentResultsDeseq_SLX14229
#Overlaid by hand for speed. For formal matching see Figure 15.
a<-    xyplot(-ma.df$log2FoldChange-0.65 ~ log(ma.df$baseMean, base=10)*0.95,
           col=c("deepskyblue4"), ylim=c(-6,2), main="Comparision between CTCF and H2av", scales="free", aspect=1, pch=20, cex=0.5,
           ylab=expression("log"[2]~"ChIP fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),
           par.settings=list(par.xlab.text=list(cex=1.1,font=2), par.ylab.text=list(cex=1.1,font=2)));
    
ma.df<-jg.experimentResultsDeseq_SLX8047

b<-    xyplot(-ma.df$log2FoldChange ~ log(ma.df$baseMean, base=10),
              col=c("palegreen3"), main="Comparision between CTCF and H2av", scales="free", aspect=1, pch=20, cex=0.5,
              ylab=expression("log"[2]~"ChIP fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),
              par.settings=list(par.xlab.text=list(cex=1.1,font=2), par.ylab.text=list(cex=1.1,font=2)));

png("plots/029_SLX-8047-SLX-14229_Overlay_CTCF_bottom.png")
a + as.layer(b)
dev.off()


png("plots/029_SLX-8047-SLX-14229_Overlay_CTCF_top.png")
b + as.layer(a)
dev.off()
