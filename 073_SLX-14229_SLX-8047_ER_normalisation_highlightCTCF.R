
######################
#
# Functions 
#
######################

setwd("/Volumes/FlyPeakCaseStudy/BrundleDevelopment")
source('package/brundle.R')

######################
#
# Settings
#
######################

dbaSummits                <- 200
jg.controlMinOverlap      <- 5
jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX14438_hs_CTCF_safe.csv"
jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX14438_hs_ER_DBA.csv"
jg.treatedCondition       =  "Fulvestrant"
jg.untreatedCondition     =  "none"

######################
#
# Main Code
#
######################

filename<-"Rdata/074_SLX-14438_dba_human_ER_CTCF.rda"
if(!file.exists(filename)){
  dbaExperiment <- jg.getDba(jg.experimentSampleSheet,dbaSummits)
  dbaControl    <- jg.getDba(jg.controlSampleSheet,   dbaSummits)
  save(dbaExperiment,dbaControl,file=filename)
} else {
  load(filename)
}

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
jg.plotMA(jg.controlPeakset,jg.experimentPeakset,jg.untreatedNames,jg.treatedNames,jg.coefficient)


#As we're going to normalise a second data set we don't need to use a DiffBind correction factor
jg.correctionFactor<-1


#Apply coefficent and control factor
jg.experimentPeaksetNormalised<-jg.applyNormalisation(jg.experimentPeakset,
                                                      jg.coefficient,
                                                      jg.correctionFactor,
                                                      jg.treatedNames
)

jg.treatedNamesPF<-jg.treatedNames
jg.untreatedNamesPF<-jg.untreatedNames

#Return values to Diffbind object
jg.dba <- DiffBind:::pv.resetCounts(dbaExperiment,
                                    jg.experimentPeaksetNormalised
)

#jg.dba is now normalised a list of ER peaks and their expected ratio


#####
#
#Rerun the normalisation code to get the ER ratio
#
######


jg.normalisedExperimentPeakset <- jg.dbaGetPeakset(jg.dba)


#Get counts for each condition
jg.ERCountsTreated<-jg.getControlCounts(jg.normalisedExperimentPeakset, 
                                             jg.experimentSampleSheet,
                                             jg.treatedCondition)
jg.ERCountsUntreated<-jg.getControlCounts(jg.normalisedExperimentPeakset,
                                          jg.experimentSampleSheet,
                                               jg.untreatedCondition)

#Get sample names for conditions
jg.treatedNames <- names(jg.ERCountsTreated)
jg.untreatedNames   <- names(jg.ERCountsUntreated)

jg.plotMA(jg.normalisedExperimentPeakset,jg.normalisedExperimentPeakset,jg.untreatedNamesPF,jg.treatedNamesPF,1)


##Establish normalised ER ratio as we know it isn't 1

jg.plotNormalization(jg.ERCountsTreated,
                     jg.ERCountsUntreated)

jg.ERratio<-jg.getNormalizationCoefficient(jg.ERCountsTreated,
                                               jg.ERCountsUntreated)

########
#
# Apply to ER only data
#
########

jg.EROnlyExperimentSampleSheet <- "samplesheet/samplesheet_SLX8047_hs.csv"

filename<-"Rdata/074_SLX-SLX-8047_dba_human_ER.rda"
if(!file.exists(filename)){
    dbaERonlyExperiment <- jg.getDba(jg.EROnlyExperimentSampleSheet,dbaSummits)
    save(dbaERonlyExperiment,file=filename)
} else {
    load(filename)
}

#Get peakset

jg.EROnlyExperimentPeakset <- jg.dbaGetPeakset(dbaERonlyExperiment)


#Get counts for each condition
jg.EROnlyCountsTreated<-jg.getControlCounts(jg.EROnlyExperimentPeakset, 
                                            jg.EROnlyExperimentSampleSheet,
                                        jg.treatedCondition)
jg.EROnlyCountsUntreated<-jg.getControlCounts(jg.EROnlyExperimentPeakset,
                                              jg.EROnlyExperimentSampleSheet,
                                          jg.untreatedCondition)

jg.plotNormalization(jg.EROnlyCountsTreated,
                     jg.EROnlyCountsUntreated)

#Get sample names for conditions
jg.treatedNames <- names(jg.EROnlyCountsTreated)
jg.untreatedNames   <- names(jg.EROnlyCountsUntreated)



jg.EROnlyExperimentRatio<-jg.getNormalizationCoefficient(jg.EROnlyCountsTreated,
                                           jg.EROnlyCountsUntreated)

#Create new ratio 
jg.ERNormalisationRatio<-jg.ERratio/jg.EROnlyExperimentRatio

#PlotMA (Optional)

jg.plotMA(jg.EROnlyExperimentPeakset,jg.EROnlyExperimentPeakset,jg.untreatedNames,jg.treatedNames,1)
jg.plotMA(jg.EROnlyExperimentPeakset,jg.EROnlyExperimentPeakset,jg.untreatedNames,jg.treatedNames,1/jg.ERNormalisationRatio)
#Value=0.8453127 

pdf("plots/pdf/074_ER-ERnormalised.pdf", points=15)
jg.plotMA(jg.EROnlyExperimentPeakset,jg.EROnlyExperimentPeakset,jg.untreatedNames,jg.treatedNames,1/jg.ERNormalisationRatio)
dev.off()


x<-jg.EROnlyExperimentPeakset
    fulvestrant<-rowMeans(x[c(4,6,8,10)])
    fulvestrant<-fulvestrant/jg.ERNormalisationRatio
    untreated<-rowMeans(x[c(5,7,9,11)])
    fc<-untreated/fulvestrant
    log2fc<- -log2(fc)
Mhs<-log2fc

Ahs<-log10(rowSums(jg.EROnlyExperimentPeakset[4:11]))

pdf("plots/pdf/074_ER-ER_normalised.pdf", points=15)
par(mar=c(5.1,5.1,4.1,2.1))
plot(Ahs,Mhs,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,pch=20,ylab=expression("log"[2]~"ChIP fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),main="Counts normalized",ylim=c(-6.25,2))
abline(h=0)
dev.off()


#NEW CODE
erChIP<-makeGRangesFromDataFrame(jg.EROnlyExperimentPeakset,keep.extra.columns=TRUE)
pfChIP<-makeGRangesFromDataFrame(jg.normalisedExperimentPeakset,keep.extra.columns=TRUE)
ctcfBED<-read.delim("SLX-14438_merged/peaks/CTCF_union.bed",header = F)
colnames(ctcfBED)<-c( "CHR" ,"START","END")
jg.controlPeaksetWide<-ctcfBED
jg.controlPeaksetWide[2]<-ctcfBED[2]-1000
jg.controlPeaksetWide[3]<-ctcfBED[3]+1000
ctcfChIP<-makeGRangesFromDataFrame(jg.controlPeaksetWide)

#jg.plotMA(jg.normalisedExperimentPeakset,jg.normalisedExperimentPeakset,jg.untreatedNamesPF,jg.treatedNamesPF,1)
#jg.plotMA(jg.EROnlyExperimentPeakset,jg.EROnlyExperimentPeakset,jg.treatedNames,jg.untreatedNames,1/jg.ERNormalisationRatio)

erChIPctcfOverlap<-subsetByOverlaps(erChIP,ctcfChIP)
pfChIPctcfOverlap<-subsetByOverlaps(pfChIP,ctcfChIP)

erChIPpfOverlap<-subsetByOverlaps(erChIPctcfOverlap,pfChIPctcfOverlap)
pfChIPerOverlap<-subsetByOverlaps(pfChIPctcfOverlap,erChIPctcfOverlap)

erChIPuntreated<-rowMeans(as.matrix(elementMetadata(erChIPpfOverlap)[,jg.untreatedNames]))
erChIPtreated<-rowMeans(as.matrix(elementMetadata(erChIPpfOverlap)[,jg.treatedNames]))
erChIPtreated/erChIPuntreated

pfChIPuntreated<-rowMeans(as.matrix(elementMetadata(pfChIPerOverlap)[,jg.untreatedNamesPF]))
pfChIPtreated<-rowMeans(as.matrix(elementMetadata(pfChIPerOverlap)[,jg.treatedNamesPF]))
pfChIPtreated/pfChIPuntreated

png("plots/073_CTCFproximal.png",pointsize=15)
boxplot(log(erChIPtreated/erChIPuntreated),
        log(pfChIPtreated/pfChIPuntreated),
        xlab="Experiment",
        ylab="Log(Fold-Change)",
        main=paste("Fold-Change of ER binding","Proximal to CTCF",collapse="\n"),
        names=c(" Cross-Normalised"," Parrallel Factor"),
        pch=20)

dev.off()
