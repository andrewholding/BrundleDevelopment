library(Brundle)

#"SLX-15439","D701-D501","ATTACTCG-TATAGCCT","PDX01-AB555B-x2" ER 8, PR 2, FOXA1 8
#"SLX-15439","D703-D502","CGCTCATT-ATAGAGGC","PDX02-VHI0244o2-x4" ER 8, PR 7, FOXA1 8
#"SLX-15439","D701-D503","ATTACTCG-CCTATCCT","PDX03-AB580-x1"  - ER 8, PR/FOXA1 unknown
#"SLX-15439","D702-D501","TCCGGAGA-TATAGCCT","PDX04-STG195-x4" - ER/PR/FOXA1/8
#"SLX-15439","D702-D503","TCCGGAGA-CCTATCCT","PDX05-V0980U-x2"  ER 5, PR 0, FOXA


pdx<-dba(sampleSheet="samplesheet/SLX-15439.csv")
pdx<-dba.count(pdx, bParallel=FALSE,bRemoveDuplicates=TRUE)
dba.plotPCA(pdx)
dba.plotHeatmap(pdx)
saveRDS(file="Rdata/pdx.rds",pdx)


ctcf<-dba(sampleSheet="samplesheet/SLX-15439_CTCF.csv")
ctcf<-dba.count(ctcf, bParallel=FALSE,minOverlap=3,bRemoveDuplicates=TRUE)
dba.plotPCA(ctcf)
dba.plotHeatmap(ctcf)
saveRDS(file="Rdata/ctcf.rds",ctcf)

pdx<-readRDS(file="Rdata/pdx.rds")
ctcf<-readRDS(file="Rdata/ctcf.rds")

library(ChIPQC)
samplesheet<-read.csv("samplesheet/SLX-15439_CTCF.csv")
exp = ChIPQC(ctcf)
ChIPQCreport(exp,reportFolder="ChIPQCreport/SLX-15439.csv")

write.table(file="txt/PDXconsensus.txt",as.data.frame(dba.peakset(pdx,bRetrieve=TRUE))[,1:3],quote=FALSE, row.names = FALSE, col.names=FALSE, sep="\t")

