

######################
#
# Functions
#
######################

setwd("/Volumes/FlyPeaks/FlyPeaks")
source('package/brundle.R')

######################
#
# Settings
#
######################

dbaSummits                <- 200
jg.controlSampleSheet     <-
    "samplesheet/samplesheet_SLX14229_hs_CTCF_DBA.csv"
jg.controlSampleSheet_1     <-
    "samplesheet/leaveOneOut/samplesheet_SLX14229_hs_CTCF_DBA_1.csv"
jg.controlSampleSheet_2     <-
    "samplesheet/leaveOneOut/samplesheet_SLX14229_hs_CTCF_DBA_2.csv"
jg.controlSampleSheet_3     <-
    "samplesheet/leaveOneOut/samplesheet_SLX14229_hs_CTCF_DBA_3.csv"
jg.experimentSampleSheet  <-
    "samplesheet/samplesheet_SLX14229_hs_ER_DBA.csv"
jg.treatedCondition       =  "Fulvestrant"
jg.untreatedCondition     =  "none"

######################
#
# Main Code
#
######################

filename <- "Rdata/051_dba.rda"
if (!file.exists(filename)) {
    jg.dba <- jg.normalize(
        jg.experimentSampleSheet,
        jg.controlSampleSheet,
        dbaSummits,
        jg.treatedCondition,
        jg.untreatedCondition
    )
    save(jg.dba, file = filename)
} else{
    load(filename)
}
#Analyze and plot with Diffbind
jg.dba_analysis <- dba.analyze(jg.dba)
#png("plots/051_SLX-14229_DiffBind_Analysis.png")
dba.plotMA(jg.dba_analysis, bFlip = TRUE)
#dev.off()

filename <- "Rdata/051_dba_1.rda"
if (!file.exists(filename)) {
    jg.dba_1 <- jg.normalize(
        jg.experimentSampleSheet,
        jg.controlSampleSheet,
        dbaSummits,
        jg.treatedCondition,
        jg.untreatedCondition
    )
    save(jg.dba_1, file = filename)
} else{
    load(filename)
}

filename <- "Rdata/051_dba_2.rda"
if (!file.exists(filename)) {
    jg.dba_2 <- jg.normalize(
        jg.experimentSampleSheet,
        jg.controlSampleSheet,
        dbaSummits,
        jg.treatedCondition,
        jg.untreatedCondition
    )
    save(jg.dba_2, file = filename)
} else{
    load(filename)
}

filename <- "Rdata/051_dba_3.rda"
if (!file.exists(filename)) {
    jg.dba_3 <- jg.normalize(
        jg.experimentSampleSheet,
        jg.controlSampleSheet,
        dbaSummits,
        jg.treatedCondition,
        jg.untreatedCondition
    )
    save(jg.dba_3, file = filename)
} else{
    load(filename)
}

system("git add . && git commit -m "auto commit" && git push")
