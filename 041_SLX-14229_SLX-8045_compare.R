

######################
#
# Functions
#
######################

setwd("/Volumes/FlyPeaks/FlyPeaks")
source('package/brundle.R')

######################
#
# Main Code
#
######################

filename <- "Rdata/041_SLX-14229_dba_human_ER_CTCF_normalised.rda"
if (!file.exists(filename)) {
    dbaSummits                <- 200
    jg.controlMinOverlap      <- 5
    jg.controlSampleSheet     <-
        "samplesheet/samplesheet_SLX14229_hs_CTCF_DBA.csv"
    jg.experimentSampleSheet  <-
        "samplesheet/samplesheet_SLX14229_hs_ER_DBA.csv"
    jg.treatedCondition       =  "Fulvestrant"
    jg.untreatedCondition     =  "none"

    filename_dba <- "Rdata/041_SLX-14229_dba_human_ER_CTCF.rda"
    if (!file.exists(filename)) {
        dbaExperiment <- jg.getDba(jg.experimentSampleSheet, dbaSummits)
        dbaControl    <- jg.getDba(jg.controlSampleSheet,   dbaSummits)
        save(dbaExperiment, dbaControl, file = filename)
    } else {
        load(filename)
    }


    #Condensed code from previous script
    jg.sampleIds <- jg.getSampleIds(jg.controlSampleSheet)
    jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
    jg.controlPeakset    <- jg.dbaGetPeakset(dbaControl)
    jg.controlCountsTreated <- jg.getControlCounts(jg.controlPeakset,
                                                   jg.controlSampleSheet,
                                                   jg.treatedCondition)
    jg.controlCountsUntreated <- jg.getControlCounts(jg.controlPeakset,
                                                     jg.controlSampleSheet,
                                                     jg.untreatedCondition)
    jg.untreatedNames <- names(jg.controlCountsUntreated)
    jg.treatedNames   <- names(jg.controlCountsTreated)
    jg.coefficient <-
        jg.getNormalizationCoefficient(jg.controlCountsTreated,
                                       jg.controlCountsUntreated)
    jg.correctionFactor <-
        jg.getCorrectionFactor(jg.experimentSampleSheet,
                               jg.treatedNames,
                               jg.untreatedNames)
    jg.experimentPeaksetNormalised <-
        jg.applyNormalisation(jg.experimentPeakset,
                              jg.coefficient,
                              jg.correctionFactor,
                              jg.treatedNames)
    jg.dba_SLX14229 <- DiffBind:::pv.resetCounts(dbaExperiment,
                                                 jg.experimentPeaksetNormalised)

    save(jg.dba_SLX14229,dbaExperiment, file = filename)
} else {
    load(filename)
}

jg.dba_analysis_SLX14229_unnormalised <- dba.analyze(dbaExperiment)
jg.dba_analysis_SLX14229 <- dba.analyze(jg.dba_SLX14229)

png("plots/041_SLX-14229_DiffBind_Analysis.png")
dba.plotMA(jg.dba_analysis_SLX14229, bFlip = TRUE)
dev.off()

png("plots/041_SLX-14229_DiffBind_Analysis_raw.png")
dba.plotMA(jg.dba_analysis_SLX14229_unnormalised, bFlip = TRUE)
dev.off()

unnormalised <-
    as.data.frame(dba.report(jg.dba_analysis_SLX14229_unnormalised, th = 1))
normalised <- as.data.frame(dba.report(jg.dba_analysis_SLX14229, th = 1))

png("plots/041_SLX-14229_p-value_comparision.png", width=480, height=480)
plot(
    unnormalised$p.value,
    normalised$p.value,
    pch = 20,
    xlab = "Unnormalised p-value",
    ylab = "Normalised p-value",
    main = "Effect of Normalisation on p-value"
)
dev.off()
