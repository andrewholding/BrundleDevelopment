setwd("/Volumes/FlyPeaks/FlyPeaks")

#install custom diffbind
#install.packages("../DiffBind_2.5.6.tar.gz", repos = NULL, type="source")

library(DiffBind)



filename<-"Rdata/024_SLX-14229_dba_human_ER.rda"
if(!file.exists(filename)){
  dba<-dba(sampleSheet = "samplesheet/samplesheet_SLX14229_hs_ER_DBA.csv")
  save(dba,file=filename)
} else {
  load(filename)
}

filename<-"Rdata/024_SLX-14229_dba.counts_human_ER.rda"
if(!file.exists(filename)){
  dba <- dba.count(dba, summits=200)
  
  dba <- dba.count(dba, peaks=NULL, score=DBA_SCORE_READS)
  
  png("plots/024_SLX-14229_dba.counts_human_ER.png")
  plot(dba)
  dev.off()
  
  save(dba,file=filename)
} else {
  load(filename)
}

dba_analysis<-dba.analyze(dba)
dba.plotMA(dba_analysis)


hsconsensus <- dba.peakset(dba_analysis, bRetrieve = T, DataType = DBA_DATA_FRAME)
#Correct names to samplesheet and test on intial data
names(hsconsensus)<-c("CHR","START","END","1a","1b","2a","2b","3b","3a")
newDBA <- DiffBind:::pv.resetCounts(dba, hsconsensus)
newDBA_analysis<-dba.analyze(newDBA)
png("plots/024_SLX-14229_MA_1_before_local_normalisation_human_ER.png")
dba.plotMA(newDBA_analysis)
dev.off()
#should give same result.


filename<-"Rdata/020_SLX-14229_hsconsensus_localnorm_ER.rda"
load(filename)


#convert back to integer +1 to avoid zeros 
hsconsensus_localnorm_ER[-c(1:3)]<-round(hsconsensus_localnorm_ER[-c(1:3)])
names(hsconsensus_localnorm_ER)<-c("CHR","START","END","1a","1b","2a","2b","3b","3a")

#rerun on local normalised data


newDBA_localnorm_ER <- DiffBind:::pv.resetCounts(dba, hsconsensus_localnorm_ER)
newDBA_analysis_localnorm_ER<-dba.analyze(newDBA_localnorm_ER)
png("plots/024_SLX-14229_MA_2_after_local_normalisation_human_ER.png")
dba.plotMA(newDBA_analysis_localnorm_ER)
dev.off()
