
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


## Changing font sizes for manuscript

jg.plotDeSeq<-function(ma.df, filename = 'file.name', p = 0.01, title.main = "Differential ChIP",log2fold =0.5, flip=FALSE)
{;
    
    if (flip == TRUE)
    {
        ma.df$log2FoldChange <- -ma.df$log2FoldChange
    }
    par(mar=c(5.1,5.1,4.1,2.1))
    xyplot(ma.df$log2FoldChange ~ log(ma.df$baseMean, base=10),#xlim=c(0,4),ylim=c(-3,1.25),
           groups=(ma.df$padj < p & abs(ma.df$log2FoldChange) > log2fold & !is.na(ma.df$padj)),
           col=c("black","red"), main=title.main, scales="free", aspect=1, pch=20, cex=0.5,
           ylab=expression("log"[2]~"ChIP fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),
           par.settings=list(axis.text=list(cex=1.5,font=1),par.main.text=list(cex=1.5,font=2),par.xlab.text=list(cex=1.5,font=2), par.ylab.text=list(cex=1.5,font=2))
           );
    
}

jg.plotDeSeqCombined <- function(jg.controlResultsDeseq,jg.experimentResultsDeseq,title.main,padjX,flip=FALSE)
{
    jg.controlResultsDeseq$group = 'a'
    jg.experimentResultsDeseq$group = 'b'
    
    if (flip == TRUE)
    {
        jg.controlResultsDeseq$log2FoldChange <- -jg.controlResultsDeseq$log2FoldChange
        jg.experimentResultsDeseq$log2FoldChange <- -jg.experimentResultsDeseq$log2FoldChange
    }
    
    for (i in 1:length(jg.experimentResultsDeseq$group)) {
        if (!is.na(jg.experimentResultsDeseq$padj[i]) & !is.na(jg.experimentResultsDeseq$log2FoldChange[i]) & jg.experimentResultsDeseq$padj[i] < padjX & jg.experimentResultsDeseq$log2FoldChange[i] < 0) {
            jg.experimentResultsDeseq$group[i] <- 'd'
        }
        else if (!is.na(jg.experimentResultsDeseq$padj[i]) & !is.na(jg.experimentResultsDeseq$log2FoldChange[i]) & jg.experimentResultsDeseq$padj[i] < padjX & jg.experimentResultsDeseq$log2FoldChange[i] > 0) {
            jg.experimentResultsDeseq$group[i] <- 'c'
        }
    }
    
    for (i in 1:length(jg.controlResultsDeseq$group)) {
        if (!is.na(jg.controlResultsDeseq$padj[i]) & !is.na(jg.controlResultsDeseq$log2FoldChange[i]) & jg.controlResultsDeseq$padj[i] < padjX & jg.controlResultsDeseq$log2FoldChange[i] < 0) {
            jg.controlResultsDeseq$group[i] <- 'f'
        }
        else if (!is.na(jg.controlResultsDeseq$padj[i]) & !is.na(jg.controlResultsDeseq$log2FoldChange[i]) & jg.controlResultsDeseq$padj[i] < padjX & jg.controlResultsDeseq$log2FoldChange[i] > 0) {
            jg.controlResultsDeseq$group[i] <- 'e'
        }
    }
    
    full.res = rbind(jg.controlResultsDeseq, jg.experimentResultsDeseq)
    par(mar=c(5.1,5.1,4.1,2.1))
    xyplot(full.res$log2FoldChange ~ log(full.res$baseMean, base=10), data = full.res,
           groups=full.res$group,
           col=c("grey40","grey80",  "#ff5454","#5480ff",  "#08298a","#750505"),
           ylab = expression('log'[2]*' Differential ChIP'),
           xlab = expression("log"[10]~"Mean of Normalized Counts"),
           aspect=1.0,
           pch=16,
           cex=0.5,
           main=title.main,
           scales=list(x=list(cex=1.5, relation = "free"), y =list(cex=1.5, relation="free")),
           between=list(y=0.5, x=0.5),
           auto.key = TRUE,
           par.settings=list(axis.text=list(cex=1.5,font=1),par.main.text=list(cex=1.5,font=2),par.xlab.text=list(cex=1.5,font=2), par.ylab.text=list(cex=1.5,font=2)),
           key=list(corner=c(1,0),
                    cex=1.0,
                    points=list(col=c( "gray80","gray40", "#ff5454", "#5480ff", "#750505", "#08298a","white"), pch=20),
                    text=list(c("Target Binding","Control Binding", "Target Decreased","Target Increased","Control Decreased","Control Increased", " "))
           ))
          
        
}

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
par(mar=c(5.1,5.1,4.1,2.1))
  jg.plotDeSeq(jg.controlResultsDeseq,
               p=0.01, 
               title.main="Fold-change in H2av binding",
               flip=T
               )
dev.off()
png("plots/029_SLX-8047_DeSeq_Analysis_ER.png")
par(mar=c(5.1,5.1,4.1,2.1))
jg.plotDeSeq(jg.experimentResultsDeseq,
             p=0.01,
             title.main="Fold-change in ER binding",
             flip=T
             )
dev.off()

#Repeat not using out control
jg.experimentDeSeqInternal<-jg.runDeSeq(jg.experimentPeaksetDeSeq, jg.conditions, NULL)
jg.experimentResultsDeseqInternal   = results(jg.experimentDeSeqInternal)

png("plots/029_SLX-8047_DeSeq_Analysis_ER_nocorrection.png")
    par(mar=c(5.1,5.1,4.1,2.1))
  jg.plotDeSeq(jg.experimentResultsDeseqInternal,
               title.main="Fold-change in ER binding (no correction)",
               p=0.01,
               flip=T
               )
dev.off()


#Draw Combined figure.
png("plots/029_SLX-8047_DeSeq_Analysis_HsDm.png")
    par(mar=c(5.1,5.1,4.1,2.1))
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
