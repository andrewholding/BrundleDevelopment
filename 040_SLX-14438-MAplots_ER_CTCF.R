
library(DiffBind)

setwd("/Volumes/FlyPeaks/resequence")


dba<-dba(sampleSheet = "samplesheet/samplesheet_SLX14229_hs_ER_DBA.csv")
dba <- dba.count(dba, bRemoveDuplicates=T)
dba_analyze<-dba.analyze(dba)
dba.plotMA(dba_analyze, contrast=1)
