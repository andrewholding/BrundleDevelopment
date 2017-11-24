
library(devtools)

install_github("andrewholding/Brundle")
library(Brundle)

setwd("/Volumes/FlyPeakCaseStudy/BrundleDevelopment")

jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX15091_CTCF.csv"
jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX15091_ER.csv"
jg.treatedCondition       =  "Estrogen"
jg.untreatedCondition     =  "none"

#####
#
# Main Code
#
#####

filename<-"Rdata/065_SLX-15091_dba_human_ER_CTCF.rda"
if(!file.exists(filename)){
    dbaExperiment <- jg.getDba(jg.experimentSampleSheet)
    dbaControl    <- jg.getDba(jg.controlSampleSheet)
    save(dbaExperiment,dbaControl,file=filename)
} else {
    load(filename)
}

jg.dba<-Brundle(
            dbaExperiment,
            dbaControl,
            jg.treatedCondition,
            jg.untreatedCondition,
            jg.experimentSampleSheet,
            jg.controlSampleSheet )
    
jg.dba_analysis<-dba.analyze(jg.dba)
dba.plotMA(jg.dba_analysis,bFlip=TRUE,th=0.05)

#Check CTCF 
ctcf.dba.analysis<-dba.analyze(dbaControl)
dba.plotMA(ctcf.dba.analysis,bFlip=TRUE) 
