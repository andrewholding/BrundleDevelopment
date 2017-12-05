library("Brundle")

jg.controlSampleSheet     <- "samplesheet/samplesheet_SLX14438_hs_CTCF_DBA.csv"
jg.experimentSampleSheet  <- "samplesheet/samplesheet_SLX14438_hs_ER_DBA.csv"
jg.treatedCondition       =  "Fulvestrant"
jg.untreatedCondition     =  "none"

setwd("/Volumes/FlyPeakCaseStudy/BrundleDevelopment")
filename<-"Rdata/048_SLX-14438_dba_human_ER_CTCF.rda"
load(filename)



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


jg.coefficient_baseline<-jg.getNormalizationCoefficient(jg.controlCountsTreated,
                                               jg.controlCountsUntreated)
rows<-nrow(jg.controlCountsTreated)

jg.coefficient<-matrix(nrow=100,ncol=100)


for (n in 100:1)
{
    for (i in 1:100)
    {
        set.seed(i)
        sampleRows<-sample(rows, rows*(n/100) )
        jg.coefficient[n,i]<-jg.getNormalizationCoefficient(jg.controlCountsTreated[sampleRows, ],
                               jg.controlCountsUntreated[sampleRows, ])
    }
    }
}

png("plots/069_stability.png",pointsize="15")
plot( rowMeans(jg.coefficient)/jg.coefficient_baseline*100-100,ylim=c(-2,+2),
      xlim=c(100,0),
      pch=20,col="blue",
      type="l",
      xlab="Percent of CTCF sites",
      ylab="Percent Error in Normalisation Coefficient",
      main="Stablity of Normalisation Coefficient",lwd=2)
lines( apply(t(jg.coefficient), 2, max) /jg.coefficient_baseline*100-100,col="black",lty=2,lwd=2) #max=1.956436
lines( apply(t(jg.coefficient), 2, min)/jg.coefficient_baseline*100-100, col="black",lty=2,lwd=2) #min=-1.970216
lines( apply(t(jg.coefficient), 2, quantile)['25%',]/jg.coefficient_baseline*100-100 ,col="black", lwd=2)#min=-0.4968219
lines( apply(t(jg.coefficient), 2, quantile)['75%',]/jg.coefficient_baseline*100-100, col="black",lwd=2)#max=-0.5254374
legend("topleft", legend=c("Average", "25%/75% quantile","Max/min error"),
       col=c("blue", "black","black"), lty=1:1:2, cex=0.8)
dev.off()
