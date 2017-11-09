setwd("/Volumes/FlyPeaks/FlyPeaks")
source('package/brundle.R')
jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX14438_hs_CTCF_DBA_noER.csv"
jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX14438_hs_ER_DBA.csv"
jg.treatedCondition       =  "Fulvestrant"
jg.untreatedCondition     =  "none"
dbaSummits                <- 200

#Fix as my computer is crashing on bParallel at the moment.
jg.getDba.linear<-function (jg.experimentSampleSheet,dbaSummits, ...)
{
    
    dba <- dba(sampleSheet = jg.experimentSampleSheet)
    if(exists("dbaSummits"))
    {
        dba <- dba.count(dba, summits=dbaSummits, ...)
    } else {
        dba <- dba.count(dba)
    }
    dba <- dba.count(dba, peaks=NULL, score=DBA_SCORE_READS)
    return(dba)
}

#Load data used in the paper.
filename<-"Rdata/040_SLX-14438_dba_human_ER_CTCF.rda"
load(filename)

#filename<-"Rdata/059_SLX-14438_dba_human_ER_CTCF_noER.rda"
#if(!file.exists(filename)){
#    dbaExperiment <- jg.getDba(jg.experimentSampleSheet,dbaSummits)
#    dbaControl    <- jg.getDba(jg.controlSampleSheet,   dbaSummits, minOverlap=5)
#    save(dbaExperiment,dbaControl,file=filename)
##} else {
#    load(filename)
#}

#dbaControl    <- jg.getDba.linear(jg.controlSampleSheet,  dbaSummits, minOverlap=5)

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
jg.plotMA(jg.experimentPeakset,jg.controlPeakset,jg.untreatedNames,jg.treatedNames,1)

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

#Run on DeSeq
jg.controlPeakset    <- jg.dbaGetPeakset(dba.count(dbaControl, peaks=NULL, score=DBA_SCORE_READS))
jg.controlPeaksetDeSeq<-jg.convertPeakset(jg.controlPeakset )
jg.controlSizeFactors = estimateSizeFactorsForMatrix(jg.controlPeaksetDeSeq)

jg.conditions <- read.csv(file=jg.controlSampleSheet, header=TRUE, sep=",")['Condition']

jg.controlDeSeq<-jg.runDeSeq(jg.controlPeaksetDeSeq, jg.conditions,jg.controlSizeFactors)
jg.controlResultsDeseq   = results(jg.controlDeSeq)

jg.plotDeSeq(jg.controlResultsDeseq,
             p=0.01, 
             title.main="Fold-change in CTCF binding",
             flip=T
)


###################
#Load in the samples with out ER peaks
#dbaControl_noER  <- jg.getDba.linear(jg.controlSampleSheet, bRemoveDuplicates=TRUE)

filename<-"Rdata/059_SLX-14438_dba_human_ER_CTCF_noER.rda"
#save(dbaControl_noER,file=filename)
load(filename)
## Repeat above analysis

jg.controlPeakset    <- jg.dbaGetPeakset(dbaControl_noER)


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
jg.plotMA(jg.experimentPeakset,jg.controlPeakset,jg.untreatedNames,jg.treatedNames,1)

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
png("plots/059_CTCF_noER_MAplot.png")
par(mar=c(5.1,5.1,4.1,4.1))
dba.plotMA(jg.dba_analysis,bFlip=TRUE,cex.lab=1.5, cex.axis=1.5, cex.main=1.25, cex.sub=1.5)
dev.off()
#Check CTCF 
ctcf.dba.analysis<-dba.analyze(dbaControl_noER)

dba.plotMA(ctcf.dba.analysis,bFlip=TRUE) 


#Run on DeSeq
jg.controlPeakset    <- jg.dbaGetPeakset(dba.count(dbaControl_noER, peaks=NULL, score=DBA_SCORE_READS))
jg.controlPeaksetDeSeq<-jg.convertPeakset(jg.controlPeakset )
jg.controlSizeFactors = estimateSizeFactorsForMatrix(jg.controlPeaksetDeSeq)

jg.conditions <- read.csv(file=jg.controlSampleSheet, header=TRUE, sep=",")['Condition']

jg.controlDeSeq<-jg.runDeSeq(jg.controlPeaksetDeSeq, jg.conditions,jg.controlSizeFactors)
jg.controlResultsDeseq   = results(jg.controlDeSeq)



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
           par.settings=list(axis.text=list(cex=1.5,font=1),par.main.text=list(cex=1.5,font=2),par.xlab.text=list(cex=1.5,font=2), par.ylab.text=list(cex=1.5,font=2)
           )
    );
    
}

png("plots/059_CTCF_noER_control_MAplot.png",pointsize=35)
jg.plotDeSeq(jg.controlResultsDeseq,
             p=0.01, 
             title.main="Fold-change in CTCF binding",
             flip=T
)
dev.off()

png("plots/059_CTCF_noER_Normalisation_plot.png",pointsize=15)
jg.plotNormalization(jg.controlCountsTreated,jg.controlCountsUntreated)
dev.off()
