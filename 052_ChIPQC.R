
setwd("/Volumes/FlyPeaks/FlyPeaks")

library(ChIPQC)

qc.samplesheets=c(
    "samplesheet/samplesheet_SLX8047.csv",
    "samplesheet/samplesheet_SLX12998.csv",
    "samplesheet/samplesheet_SLX14229.csv",
    "samplesheet/samplesheet_SLX14438_hs_DBA.csv",
    "samplesheet/samplesheet_SLX15090.csv",
    "samplesheet/samplesheet_SLX15091.csv")

#library(BiocParallel)
#register(MulticoreParam())
#register(SerialParam())

for (qc.samplesheet in qc.samplesheets) {
    dba.ChIPQC<-ChIPQC(qc.samplesheet, "hg19", bCount=T)
    ChIPQCreport(dba.ChIPQC,facetBy="Condition",facet=T,
                 reportFolder=paste0("ChIPQCreport/", substr(qc.samplesheet,13,100))
                )
}
