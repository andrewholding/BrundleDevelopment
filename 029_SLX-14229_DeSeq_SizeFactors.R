
library(DiffBind)
library(Rsamtools)
library(codetools)
library(DESeq2)
library(lattice)

######################
#
# Functions 
#
######################

jg.countAlignedMReads<- function(jg.bamFiles){
  jg.counts<-numeric()
  for (jg.bam in jg.bamFiles){
    jg.counts<-append(jg.counts,colSums(idxstatsBam(jg.bam)["mapped"])/1E6)
  }
  return(jg.counts)
}


jg.getControlCounts <- function(jg.control,jg.controlSampleSheet,jg.Condition)
{
  jg.controlCounts<-jg.control[,-c(1:3)]
  temp <- read.csv(file=jg.controlSampleSheet, header=TRUE, sep=",")['Condition']==jg.Condition
  return(jg.controlCounts[,temp])
}

jg.plotNormalization<-function(jg.controlCountsTreated,jg.controlCountsUntreated)
{
  plot(rowMeans(jg.controlCountsTreated),rowMeans(jg.controlCountsUntreated), pch=20,
       xlab="Counts in peak after treatment" ,  ylab="Counts in peak before treatment" ,
       main="Comparision of Counts in peaks for CTCF")
  lm1<-lm(rowMeans(jg.controlCountsUntreated) ~ 0 + rowMeans(jg.controlCountsTreated))

  abline(c(0,lm1$coef),col="red3")
  print(lm1$coefficients)
  angularcoeff<-lm1$coef[1]
  
  points(rowMeans(jg.controlCountsTreated)*angularcoeff,rowMeans(jg.controlCountsUntreated),pch=20, col="royalblue3" )
  treatment_fit<-rowMeans(jg.controlCountsTreated)*angularcoeff
  lm1<-lm(treatment_fit ~ 0 + rowMeans(jg.controlCountsUntreated))
  abline(c(0,lm1$coef),col="purple")
  legend("topleft",legend=c("Raw", "Normalised"),pch=20,col=c("black","royalblue3"))
}


jg.getNormalizationCoefficient<-function(jg.controlCountsTreated,jg.controlCountsUntreated)
{
  lm1<-lm(rowMeans(jg.controlCountsUntreated) ~ 0 + rowMeans(jg.controlCountsTreated))
  return(lm1$coef[1])
}

jg.plotMA<-function(jg.experimentPeakset,jg.controlPeakset,jg.untreatedNames,jg.treatedNames,jg.coefficient)
{
  M_corrected<-apply(jg.experimentPeakset[-c(1:3)],1,function(x){
    untreated<-mean(x[jg.untreatedNames])
    treated<-jg.coefficient*mean(x[jg.treatedNames])
    fc<-mean(treated)/mean(untreated)
    log2fc<-log2(fc)
    return(log2fc)
  })
  A_corrected<-apply(jg.experimentPeakset[-c(1:3)],1,function(x){
    untreated<-mean(x[jg.untreatedNames])
    treated<-jg.coefficient*mean(x[jg.treatedNames])
    return(log10(sum(treated+untreated)))
  })
  
  
  M_dm_corrected<-apply(jg.controlPeakset[-c(1:3)],1,function(x){
    untreated<-mean(x[jg.untreatedNames])
    treated<-jg.coefficient*mean(x[jg.treatedNames])
    fc<-mean(treated)/mean(untreated)
    log2fc<-log2(fc)
    return(log2fc)
  })
  
  
  A_dm_corrected<-apply(jg.controlPeakset[-c(1:3)],1,function(x){
    untreated<-mean(x[jg.untreatedNames])
    treated<-jg.coefficient*mean(x[jg.treatedNames])
    return(log10(sum(treated+untreated)))
  })
  
  plot(A_corrected,M_corrected,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(treatment)", main="Normalised aligned reads")
  points(A_dm_corrected,M_dm_corrected,pch=20,col="cornflowerblue")
  lm1<-lm(M_dm_corrected~A_dm_corrected)
  abline(lm1$coef,col="red4")
  abline(h=0)
  
}

jg.getDba<-function (jg.experimentSampleSheet,dbaSummits)
{
  
  dba <- dba(sampleSheet = jg.experimentSampleSheet)
  dba <- dba.count(dba, summits=dbaSummits)
  dba <- dba.count(dba, peaks=NULL, score=DBA_SCORE_READS)
  return(dba)
}


jg.dbaGetPeakset <-function(dba)
{
  jg.peakset<-dba.peakset(dba, bRetrieve = T, DataType = DBA_DATA_FRAME)
  #Correct sample names back to that in sample sheet, 
  #as DiffBind changes them on export.

  jg.sampleIds<-dba$samples[,'SampleID']
  names(jg.peakset)<-c("CHR","START","END",jg.sampleIds)
  return(jg.peakset)
}

jg.getSampleIds<-function(jg.controlSampleSheet)
{
  jg.sampleIds<-as.character(read.csv(file=jg.controlSampleSheet, header=TRUE, sep=",")[,'SampleID'])
  return(jg.sampleIds)
}

jg.getCorrectionFactor <-function (jg.experimentSampleSheet,jg.treatedNames,jg.untreatedNames)
{
  #Load aligned reads for experiment Bams. 
  
  #Get list of experimentBams
  jg.experimentBams <- as.character(read.csv(file=jg.experimentSampleSheet, header=TRUE, sep=",")[,'bamReads'])
  jg.sampleIds<-as.character(read.csv(file=jg.experimentSampleSheet, header=TRUE, sep=",")[,'SampleID'])
  jg.experimentAligned<-jg.countAlignedMReads(jg.experimentBams)
  names(jg.experimentAligned)<-jg.sampleIds
  #Take ratio of treated:untreated aligned reads and create correction factor.
  jg.correctionFactor<-sum(jg.experimentAligned[jg.treatedNames])/sum(jg.experimentAligned[jg.untreatedNames])
  return(jg.correctionFactor) 
}

jd.applyNormalisation<-function(jg.experimentPeakset,jg.coefficient, jg.correctionFactor,jg.treatedNames)
{
  jg.experimentPeaksetNormalised<-jg.experimentPeakset
  jg.experimentPeaksetNormalised[jg.treatedNames]<-
    (jg.coefficient*jg.correctionFactor*jg.experimentPeakset[jg.treatedNames])
  return(jg.experimentPeaksetNormalised)
}


plot.ma.lattice <- function(ma.df, filename = 'file.name', p = 0.01, title.main = "Differential ChIP",log2fold =0.5)
{;
  #Adapted from ....
  #pdf(paste("MA_plot", filename, ".pdf", sep=''), width=3.83, height=3.83);
  print(xyplot(ma.df$log2FoldChange ~ log(ma.df$baseMean, base=10),
               groups=(ma.df$padj < p & abs(ma.df$log2FoldChange) > log2fold & !is.na(ma.df$padj)),
               col=c("black","red"), main=title.main, scales="free", aspect=1, pch=20, cex=0.5,
               ylab=expression("log"[2]~"ChIP fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),
               par.settings=list(par.xlab.text=list(cex=1.1,font=2), par.ylab.text=list(cex=1.1,font=2))));
  #dev.off()
}

#Check Functions
checkUsage(jg.countAlignedMReads)
checkUsage(jg.getControlCounts)
checkUsage(jg.plotNormalization)
checkUsage(jg.getNormalizationCoefficient)
checkUsage(jg.MAplot) 
checkUsage(jg.getDba)
checkUsage(jg.dbaGetPeakset)
checkUsage(jg.getSampleIds)
checkUsage(jg.getCorrectionFactor)
checkUsage(jd.applyNormalisation)
######################
#
# Settings
#
######################

setwd("/Volumes/FlyPeaks/FlyPeaks")
dbaSummits                <- 200
jg.controlMinOverlap      <- 5
jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX14229_hs_CTCF_DBA.csv"
jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX14229_hs_ER_DBA.csv"
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

#Convert Peakset to DeSeq Workflow
jg.controlPeaksetDeSeq<-jg.controlPeakset[-c(1:3)]
row.names(jg.controlPeaksetDeSeq)<-
   paste(jg.controlPeakset[,1], ':', jg.controlPeakset[,2], '-', jg.controlPeakset[,3], sep='')
  
#Establish size factors directly from Control data
jf.controlSizeFactors = estimateSizeFactorsForMatrix(jg.controlPeaksetDeSeq)

#Get conditions dataframe for DeSeq
jf.conditions <- read.csv(file=jg.controlSampleSheet, header=TRUE, sep=",")['Condition']

#Run DeSeq on control
jg.controlDeSeq = DESeqDataSetFromMatrix(jg.controlPeaksetDeSeq, jf.conditions, ~Condition)
jg.controlDeSeq = estimateSizeFactors(jg.controlDeSeq)
jg.controlDeSeq = estimateDispersions(jg.controlDeSeq)
jg.controlDeSeq = nbinomWaldTest(jg.controlDeSeq)
jg.controlResultsDeseq   = results(jg.controlDeSeq)

#Repeat for experiment

#Convert experiment Peakset to DeSeq Workflow
jg.experimentPeaksetDeSeq<-jg.experimentPeakset[-c(1:3)]
row.names(jg.experimentPeaksetDeSeq)<-
  paste(jg.experimentPeakset[,1], ':', jg.experimentPeakset[,2], '-', jg.experimentPeakset[,3], sep='')

jg.experimentDeSeq = DESeqDataSetFromMatrix(jg.experimentPeaksetDeSeq, jf.conditions, ~Condition)
sizeFactors(jg.experimentDeSeq) = jf.controlSizeFactors #Use control size factors
jg.experimentDeSeq         = estimateDispersions(jg.experimentDeSeq)
jg.experimentDeSeq = nbinomWaldTest(jg.experimentDeSeq)
jg.experimentResultsDeseq   = results(jg.experimentDeSeq)

plot.ma.lattice(jg.controlResultsDeseq, p=0.01, filename = 'test.ctcf')
plot.ma.lattice(jg.experimentResultsDeseq, p=0.01 ,filename = 'er.ctcf')

#Repeat not using out control
jg.experimentDeSeqInternal = estimateSizeFactors(jg.experimentDeSeq)
jg.experimentDeSeqInternal = nbinomWaldTest(jg.experimentDeSeqInternal)
jg.experimentResultsDeseqInternal   = results(jg.experimentDeSeqInternal)

plot.ma.lattice(jg.experimentResultsDeseqInternal, p=0.01, filename = 'test.ctcf')

#Just drawing figure
jg.controlResultsDeseq$group = 'a'
jg.experimentResultsDeseq$group = 'b'

padjX = 0.01

for (i in 1:length(jg.experimentResultsDeseq$group)) {
  if (!is.na(jg.experimentResultsDeseq$padj[i]) & !is.na(jg.experimentResultsDeseq$log2FoldChange[i]) & jg.experimentResultsDeseq$padj[i] < padjX & jg.experimentResultsDeseq$log2FoldChange[i] < 0) {
    jg.experimentResultsDeseq$group[i] <- 'c'
  }
  else if (!is.na(jg.experimentResultsDeseq$padj[i]) & !is.na(jg.experimentResultsDeseq$log2FoldChange[i]) & jg.experimentResultsDeseq$padj[i] < padjX & jg.experimentResultsDeseq$log2FoldChange[i] > 0) {
    jg.experimentResultsDeseq$group[i] <- 'd'
  }
}

for (i in 1:length(jg.controlResultsDeseq$group)) {
  if (!is.na(jg.controlResultsDeseq$padj[i]) & !is.na(jg.controlResultsDeseq$log2FoldChange[i]) & jg.controlResultsDeseq$padj[i] < padjX & jg.controlResultsDeseq$log2FoldChange[i] < 0) {
    jg.controlResultsDeseq$group[i] <- 'e'
  }
  else if (!is.na(jg.controlResultsDeseq$padj[i]) & !is.na(jg.controlResultsDeseq$log2FoldChange[i]) & jg.controlResultsDeseq$padj[i] < padjX & jg.controlResultsDeseq$log2FoldChange[i] > 0) {
    jg.controlResultsDeseq$group[i] <- 'f'
  }
}




full.res = rbind(jg.controlResultsDeseq, jg.experimentResultsDeseq)

#pdf('Holding_plot.pdf', width=3.5, height=3.5, useDingbats=FALSE)
xyplot(full.res$log2FoldChange ~ log(full.res$baseMean, base=10), data = full.res,
       groups=full.res$group,
       col=c("grey40","grey80",  "#ff5454", "#5480ff", "#750505", "#08298a"),
       ylab = expression('log'[2]*' Differential ChIP'),
       xlab = expression("log"[10]~"Mean of Normalized Counts"),
       aspect=1.0,
       pch=16,
       cex=0.5,
       scales=list(x=list(cex=0.8, relation = "free"), y =list(cex=0.8, relation="free")),
       between=list(y=0.5, x=0.5),
)
#dev.off()


