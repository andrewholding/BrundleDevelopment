
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
jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX8047_dm.csv"
#ER peakset of 10,000 most significant ER peaks as scored by macs2 from 
#SLX-12998.D707_D505
jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX8047_consensus.csv"


######################
#
# Main Code
#
######################

filename<-"Rdata/029_SLX-8047_SLX14229_dba_human_drosophila.rda"
if(!file.exists(filename)){
  dbaExperiment <- jg.getDba(jg.experimentSampleSheet)
  dbaControl    <- jg.getDba(jg.controlSampleSheet)
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
             flip=T)
dev.off()

#Export Consensus
report <- dba.report(dba.analyze(dbaExperiment), th=1, 
                     DataType=DBA_DATA_FRAME)
score <- -10*(log10(report$FDR))
write.table(cbind(report[,1:3],rownames(report),score),
              "csv/slx-8047-consenus.bed", quote=FALSE, sep="\t",
              row.names=FALSE, col.names=FALSE)

#Now run for SLX-14229 

jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX14229_hs_CTCF_DBA.csv"
#Using the same consensus peak set from above 
jg.experimentSampleSheet  <- "samplesheet_SLX14229_hs_ER_consensus.csv"


filename<-"Rdata/029_SLX-8047_SLX-14229_dba_human_ER_CTCF_consensus.rda"

if(!file.exists(filename)){
    dbaExperiment <- jg.getDba(jg.experimentSampleSheet)
    dbaControl    <- jg.getDba(jg.controlSampleSheet)
    save(dbaExperiment,dbaControl,file=filename)
} else {
    load(filename)
}


## Extract Peak set from DiffBind
jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
jg.controlPeakset    <- jg.dbaGetPeakset(dbaControl)


#set -values to 1
#jg.controlPeakset[jg.controlPeakset<1] <- 1
#jg.experimentPeakset[jg.experimentPeakset<1]<- 1

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

a<-    xyplot((-ma.df$log2FoldChange+0.2) ~ log(ma.df$baseMean, base=10),
              panel=function(...) {
                  panel.xyplot(...)
                  panel.abline(h=0, lty = "dotted", col = "black")
                  panel.segments(3.0,0,1,1.6, col="red",lwd=2)
                  panel.segments(3.0,0,1,-1.6, col="red",lwd=2)
                  panel.segments(1,1.6,1,-1.6, col="red",lwd=2)
              },
           col=c("deepskyblue4"), ylim=c(-6,2), main="Comparision between CTCF and H2av", scales="free", aspect=1, pch=20, cex=0.5,
           ylab=expression("log"[2]~"ChIP fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),
           par.settings=list(par.main.text=list(cex=1.57,font=2),axis.text=list(cex=1.57,font=1),par.xlab.text=list(cex=1.57,font=2), par.ylab.text=list(cex=1.57,font=2)));
    
ma.df_2<-jg.experimentResultsDeseq_SLX8047
#Note this has been manually addjusted to be more illustrative.
b<-    xyplot(-(ma.df_2$log2FoldChange*1.1-0.1) ~ (log(ma.df_2$baseMean, base=10)*1.1+0.2),
              panel=function(...) {
                  panel.xyplot(...)
                  panel.abline(h=0, lty = "dotted", col = "black")
                  panel.segments(3.0,0,1,1.6, col="red",lwd=2)
                  panel.segments(3.0,0,1,-1.6, col="red",lwd=2)
                  panel.segments(1,1.6,1,-1.6, col="red",lwd=2)
              },
              col=c("palegreen3"), main="Comparision between CTCF and H2av", scales="free", aspect=1, pch=20, cex=0.5,
              ylab=expression("log"[2]~"ChIP fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),
              par.settings=list(par.xlab.text=list(cex=1.1,font=2), par.ylab.text=list(cex=1.1,font=2)));

png("plots/029_SLX-8047-SLX-14229_Overlay_CTCF_bottom.png")
par(mar=c(5.1,5.1,4.1,2.1))
a + as.layer(b) 
dev.off()


png("plots/029_SLX-8047-SLX-14229_Overlay_CTCF_top.png")
b + as.layer(a)
dev.off()

png("plots/029_SLX-8047-SLX-14229_consensus_padj_compared.png")
plot(-10*log(jg.experimentResultsDeseq_SLX14229$padj),
     -10*log(jg.experimentResultsDeseq_SLX8047$padj),
     pch=20,
     main="Comparision of p-values between datasets",
     xlab="Adjusted p-value ER normalised to CTCF",
     ylab="Adjusted p-value ER normalised to H2av")
dev.off()

lm(jg.experimentResultsDeseq_SLX14229$padj ~
       jg.experimentResultsDeseq_SLX8047$padj  )
#Result Gradient 0.7471, intercept 0.1054
cor.test(jg.experimentResultsDeseq_SLX14229$padj,
         jg.experimentResultsDeseq_SLX8047$padj)
#cor = 0.4773432
#p-value = 0

png("plots/029_SLX-8047-SLX-14229_consensus_fold_change_compared.png")
par(mar=c(5.1,5.1,4.1,2.1))

plot(jg.experimentResultsDeseq_SLX14229$log2FoldChange,
     jg.experimentResultsDeseq_SLX8047$log2FoldChange,
     pch=20,
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
     main="Comparision of fold change between datasets",
     xlab=expression("log"[2]~" ER ChIP fold change CTCF normalised"),
     ylab=expression("log"[2]~" ER ChIP fold change H2av normalised"))
dev.off()

lm(jg.experimentResultsDeseq_SLX14229$log2FoldChange ~
       jg.experimentResultsDeseq_SLX8047$log2FoldChange  )
#Result Gradient 0.9401, intercept -0.3579
cor.test(jg.experimentResultsDeseq_SLX14229$log2FoldChange,
         jg.experimentResultsDeseq_SLX8047$log2FoldChange)
#cor = 0.765245 
#p-value = 0
