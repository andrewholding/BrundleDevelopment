
library(DiffBind)
library(Rsamtools)
library(codetools) 

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
       main="Comparision of Counts in peaks for Drosophila")
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

jg.MAplot<-function(jg.experimentPeakset,jg.controlPeakset,jg.untreatedNames,jg.treatedNames,jg.coefficient)
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

#Try setting Reads to unity instead
jg.correctionFactor<-1

#Apply coefficent and control factor
jg.experimentPeaksetNormalised<-jd.applyNormalisation(jg.experimentPeakset,
                                                      jg.coefficient,
                                                      jg.correctionFactor
                                                      jg.treatedNames
                                                      )
 

#Return values to Diffbind and plot normalised result.
jg.dba <- DiffBind:::pv.resetCounts(dbaExperiment,
                                    jg.experimentPeaksetNormalised
                                    )
jg.dba$class["Reads",]<-as.numeric(jg.dba$class["Reads",])*0+1  
           
#Analyze and plot with Diffbind                                                           
jg.dba_analysis<-dba.analyze(jg.dba)
dba.plotMA(jg.dba_analysis)

#Original vs Normalised data for comparison
dba_analysis<-dba.analyze(dbaExperiment)
par(mfrow=c(1,2))
dba.plotMA(dba_analysis)
dba.plotMA(jg.dba_analysis)
par(mfrow=c(1,1))

