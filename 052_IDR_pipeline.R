setwd("~/Dropbox (Cambridge University)/idr in R")

#http://web.cs.ucdavis.edu/~filkov/software/iplant/idr/functions-all-clayton-12-13.r
source("functions-all-clayton-12-13.r")
chr.size <- read.table("genome_table.txt")
##########
#
# Functions
#
##########


IDR <- function(peakfile1, peakfile2,sampleName1,sampleName2,experiment, condition)
{
    #Code from Encode
    half.width <- NULL
    overlap.ratio <- 0
    is.broadpeak <- F
    sig.value <- "p.value"
    
    rep1 <- process.narrowpeak(
        paste(peakfile1, sep = ""),
        chr.size,
        half.width = half.width,
        summit = "offset",
        broadpeak = is.broadpeak
    )
    rep2 <- process.narrowpeak(
        paste(peakfile2, sep = ""),
        chr.size,
        half.width = half.width,
        summit = "offset",
        broadpeak = is.broadpeak
    )
    uri.output <-
        compute.pair.uri(
            rep1$data.cleaned,
            rep2$data.cleaned,
            sig.value1 = sig.value,
            sig.value2 = sig.value,
            overlap.ratio = overlap.ratio
        )
    em.output <- fit.em(uri.output$data12.enrich, fix.rho2 = T)
    idr.local <- 1 - em.output$em.fit$e.z
    IDR <- c()
    o <- order(idr.local)
    IDR[o] <- cumsum(idr.local[o]) / c(1:length(o))
    idr_output <-
        data.frame(
            chr1 = em.output$data.pruned$sample1[, "chr"],
            start1 = em.output$data.pruned$sample1[, "start.ori"],
            stop1 = em.output$data.pruned$sample1[, "stop.ori"],
            sig.value1 = em.output$data.pruned$sample1[, "sig.value"],
            chr2 = em.output$data.pruned$sample2[, "chr"],
            start2 = em.output$data.pruned$sample2[, "start.ori"],
            stop2 = em.output$data.pruned$sample2[, "stop.ori"],
            sig.value2 = em.output$data.pruned$sample2[, "sig.value"],
            idr.local = 1 - em.output$em.fit$e.z,
            IDR = IDR
        )
    
    write.table(idr_output,
                "idr_overlapped_peaks.txt",
                sep = "",
                quote = F)
    
    
    filtered_peaks <- idr_output[idr_output[, 10] <= 0.01, ]
    dim(filtered_peaks) # get the number of peaks
    
    
    ###Plotting
    
    ez.list <-
        get.ez.tt.all(em.output,
                      uri.output$data12.enrich$merge1,
                      uri.output$data12.enrich$merge2)
    par(
        mar = c(5, 5, 0, 0.5),
        mfrow = c(1, 3),
        oma = c(5, 0, 2, 0)
    )
    idr_output$col[idr_output[, 10] <= 0.01] = "black"
    idr_output$col[idr_output[, 10] >= 0.01] = "red"
    plot(
        log(idr_output[, 4]),
        log(idr_output[, 8]),
        col = idr_output[, 11],
        pch = 19,
        xlab = paste0("log(signal) ", sampleName1),
        ylab = paste0("log(signal) ", sampleName2),
    )
    legend(
        "topleft",
        c("IDR=>0.01", "IDR<=0.01"),
        col = c("red", "black"),
        pch = 19,
        bty = "n",
        lty = c(1, 1),
        lwd = c(2, 2)
    )
    plot(
        rank(-idr_output[, 4]),
        rank(-idr_output[, 8]),
        col = idr_output[, 11],
        pch = 19,
        xlab = paste0("log(signal) ", sampleName1),
        ylab = paste0("log(signal) ", sampleName2),
    )
    legend(
        "topleft",
        c("IDR=>0.01", "IDR<=0.01"),
        col = c("red", "black"),
        pch = 19,
        bty = "n",
        lty = c(1, 1),
        lwd = c(1, 1)
    )
    plot(ez.list$IDR, ylab = "IDR", xlab = "num of significant peaks")
    mtext(paste(experiment," ", condition) ,side = 3, line = -51, outer = TRUE)
    
}

##########
#
# Analysis
#
#########

setwd("/Volumes/FlyPeaks/FlyPeaks")

pdf("plots/052_idr.pdf")
idr.samplesheets=c(
    "samplesheet/samplesheet_SLX8047.csv",
    "samplesheet/samplesheet_SLX12998.csv",
    "samplesheet/samplesheet_SLX14229.csv",
    "samplesheet_SLX14438_hs_DBA.csv")

for (idr.samplesheet in idr.samplesheets) {

idr.sampletable<-read.csv( idr.samplesheet)

idr.conditions<-levels(idr.sampletable$Condition)

for (condition in idr.conditions) {
    idr.samples<- idr.sampletable[idr.sampletable$Condition == condition,]$Peaks
    idr.noOfSamples<-length(idr.samples)
    idr.combinations<-combn(idr.noOfSamples,2, simplify=TRUE)
    for (combination in seq(1,length(idr.combinations)/2)) {
        idr.samples[idr.combinations[,combination]]
        sampleName1<-idr.sampletable[idr.sampletable$Peaks==idr.samples[idr.combinations[,combination]][1],]$SampleID
        sampleName2<-idr.sampletable[idr.sampletable$Peaks==idr.samples[idr.combinations[,combination]][2],]$SampleID
        IDR(idr.samples[idr.combinations[,combination]][1],idr.samples[idr.combinations[,combination]][2],sampleName1,sampleName2,idr.samplesheet,condition)
    }
}

}
dev.off()
