setwd("/Volumes/FlyPeaks/FlyPeaks")
library(DiffBind)


filename<-"Rdata/025_SLX-14229_dba_human_ER.rda"
if(!file.exists(filename)){
  dba<-dba(sampleSheet = "samplesheet/samplesheet_SLX14229_hs_ER_DBA.csv")
  save(dba,file=filename)
} else {
  load(filename)
}


filename<-"Rdata/025_SLX-14229_dba.counts_human_ER.rda"
if(!file.exists(filename)){
  dba <- dba.count(dba, summits=200)
  dba <- dba.count(dba, peaks=NULL, score=DBA_SCORE_READS)
  save(dba,file=filename)
} else {
  load(filename)
}

dba_analysis<-dba.analyze(dba)
png("plots/025_diffbindMA_import.png")
dba.plotMA(dba_analysis)
dev.off()

hsconsensus <- dba.peakset(dba_analysis, bRetrieve = T, DataType = DBA_DATA_FRAME)
#Correct names to samplesheet and test on intial data
names(hsconsensus)<-c("CHR","START","END","1a","1b","2a","2b","3b","3a")
newDBA <- DiffBind:::pv.resetCounts(dba, hsconsensus)
newDBA_analysis<-dba.analyze(newDBA)

png("plots/025_diffbindMA_reimported.png")
dba.plotMA(newDBA_analysis)
dev.off()

#Now lets try a linear model.
hsconsensus_ER<-readRDS("Rdata/018_SLX-14229_hsconsensus_ER.rds")
hsconsensus_CTCF<-readRDS("Rdata/018_SLX-14229_hsconsensus_CTCF.rds")
hscounts_ER<-hsconsensus_ER[,-c(1:3)]
hscounts_CTCF<-hsconsensus_CTCF[,-c(1:3)]


load("Rdata/012_SLX-14229_aligned.rda")
aligned<-aligned[grep("SLX-14229.D",names(aligned))]
aligned<-sapply(aligned,sum)/1E6
#remove inputs/controls to match samples here
aligned<-aligned[c(1,2,5,7,8,10)]

hsrpm_ER<-hscounts_ER
for(i in 1:length(hsrpm_ER)){
  hsrpm_ER[i]<-hscounts_ER[i]/aligned[i]
}

hsrpm_CTCF<-hscounts_CTCF
for(i in 1:length(hsrpm_CTCF)){
  hsrpm_CTCF[i]<-hscounts_CTCF[i]/aligned[i]
}



M_RPM_ER<-apply(hsrpm_ER,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])
  untreated<-mean(x[c(1,3,6)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A_RPM_ER<-apply(hsrpm_ER,1,function(x){
  return(log10(sum(x)))
})

M_RPM_CTCF<-apply(hsrpm_CTCF,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])
  untreated<-mean(x[c(1,3,6)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A_RPM_CTCF<-apply(hsrpm_CTCF,1,function(x){
  return(log10(sum(x)))
})

png("plots/025_MAplot_showing_CTCF_reads.png")
plot(A_RPM_ER,M_RPM_ER,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="RPM aligned reads")
points(A_RPM_CTCF,M_RPM_CTCF,pch=20,col="gray")
lm1<-lm(M_RPM_CTCF~A_RPM_CTCF)
abline(lm1$coef,col="blue")
abline(h=0)
dev.off()

lm1<-lm(M_RPM_CTCF~A_RPM_CTCF)
intercept<-lm1$coef[1]
angularcoeff<-lm1$coef[2]
MFit_RPM_ER<-M_RPM_ER-(A_RPM_ER*angularcoeff)-intercept
MFit_RPM_CTCF<-M_RPM_CTCF-(A_RPM_CTCF*angularcoeff)-intercept

png("plots/025_MAplot_showing_CTCF_reads.png")
plot(A_RPM_ER,MFit_RPM_ER,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="RPM aligned reads - CTCF fit linear fit to Log2 ratio")
points(A_RPM_CTCF,MFit_RPM_CTCF,pch=20,col="gray")
lm1<-lm(MFit_RPM_CTCF~A_RPM_CTCF)
abline(lm1$coef,col="blue")
abline(h=0)
dev.off()

#Non-trivial to use linear fit of ratios, there fit data before generating ratios
fulvestrant_CTCF<-rowMeans(hsrpm_CTCF[c(2,4,5)])
untreated_CTCF<-rowMeans(hsrpm_CTCF[c(1,3,6)])


png("plots/025_CTCF_peak_count_comparision_linear.png")
plot(fulvestrant_CTCF,untreated_CTCF, pch=20, xlab="Counts in peak after treatment" ,  ylab="Counts in peak before treatment" , main="Comparision of counts in peaks for CTCF")
lm1<-lm(untreated_CTCF~ 0 + fulvestrant_CTCF)
abline(c(0,lm1$coef),col="red3")
intercept<-0
angularcoeff<-lm1$coef[1] #is > 1 which implies less fuvestrant cells
fulvestrant_CTCF_fit<-fulvestrant_CTCF*angularcoeff+intercept
points(fulvestrant_CTCF_fit,untreated_CTCF,pch=20, col="royalblue3" )
legend("topleft",legend=c("Raw", "Normalised"),pch=20,col=c("black","royalblue3"))
abline(c(0,1),col="purple3")
dev.off()
#test on MA plots



M_RPM_CTCF_corrected<-apply(hsrpm_CTCF,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])*angularcoeff+intercept
  untreated<-mean(x[c(1,3,6)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A_RPM_CTCF_corrected<-apply(hsrpm_CTCF,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])*angularcoeff+intercept
  untreated<-mean(x[c(1,3,6)])
  return(log10(sum(fulvestrant+untreated)))
})



M_RPM_ER_corrected<-apply(hsrpm_ER,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])*angularcoeff
  untreated<-mean(x[c(1,3,6)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A_RPM_ER_corrected<-apply(hsrpm_ER,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])*angularcoeff
  untreated<-mean(x[c(1,3,6)])
  return(log10(sum(fulvestrant+untreated)))
})

png("plots/025_MAplot_linear_CTCF_fit.png")
plot(A_RPM_ER_corrected,M_RPM_ER_corrected,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="RPM aligned reads - linear CTCF Fit")
points(A_RPM_CTCF_corrected,M_RPM_CTCF_corrected,pch=20,col="gray")
lm1<-lm(M_RPM_CTCF_corrected~A_RPM_CTCF_corrected)
abline(lm1$coef,col="blue")
abline(h=0)
dev.off()



#Lets try correcting the log2 values

png("plots/025_CTCF_peak_count_comparision_log.png")
plot(log2(fulvestrant_CTCF),log2(untreated_CTCF), pch=20,
     xlab="Log2(counts) in peak after treatment" ,  ylab="Log2(counts) in peak before treatment" ,
     main="Comparision of Log2(counts) in peaks for CTCF")

lm1<-lm(log2(untreated_CTCF)~ 0 + log2(fulvestrant_CTCF))
abline(c(0,lm1$coef),col="red3")
abline(h=0)

intercept<-0
angularcoeff<-lm1$coef[1] #is > 1 which implies less fuvestrant cells

fulvestrant_CTCF_log2fit<-fulvestrant_CTCF**angularcoeff
points(log2(fulvestrant_CTCF_log2fit),log2(untreated_CTCF),pch=20, col="royalblue3" )
abline(c(0,lm1$coef),col="red3")
abline(h=0)
lm1<-lm(log2(untreated_CTCF)~ 0 + log2(fulvestrant_CTCF_log2fit))
abline(c(0,lm1$coef),col="purple")
abline(h=0)
legend("topleft",legend=c("Raw", "Normalised"),pch=20,col=c("black","royalblue3"))
dev.off()

log2fit_coeff<-angularcoeff

M_RPM_CTCF_log_corrected<-apply(hsrpm_CTCF,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])**angularcoeff
  untreated<-mean(x[c(1,3,6)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A_RPM_CTCF_log_corrected<-apply(hsrpm_CTCF,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])**angularcoeff
  untreated<-mean(x[c(1,3,6)])
  return(log10(sum(fulvestrant+untreated)))
})


M_RPM_ER_log_corrected<-apply(hsrpm_ER,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])**angularcoeff
  untreated<-mean(x[c(1,3,6)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A_RPM_ER_log_corrected<-apply(hsrpm_ER,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])**angularcoeff
  untreated<-mean(x[c(1,3,6)])
  return(log10(sum(fulvestrant+untreated)))
})

png("plots/025_MAplot_linear_CTCF_log2fit.png")
plot(A_RPM_ER_log_corrected,M_RPM_ER_log_corrected,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="RPM aligned reads - log2 CTCF Fit")
points(A_RPM_CTCF_log_corrected,M_RPM_CTCF_log_corrected,pch=20,col="gray")
lm1<-lm(M_RPM_CTCF_log_corrected~A_RPM_CTCF_log_corrected)
abline(lm1$coef,col="blue")
abline(h=0)
dev.off()

#Apply the log2 CTCF fit
hsconsensus <- dba.peakset(dba_analysis, bRetrieve = T, DataType = DBA_DATA_FRAME)
#Correct names to samplesheet and test on intial data
names(hsconsensus)<-c("CHR","START","END","1a","1b","2a","2b","3b","3a")
hsconsensus_CTCFfit<-hsconsensus
#change to RPM, fit and then convert back. 
hsconsensus_CTCFfit[c(5,7,8)]<-round(((hsconsensus[c(5,7,8)]/aligned[c(2,4,5)])**log2fit_coeff)*aligned[c(2,4,5)])
  
newDBA <- DiffBind:::pv.resetCounts(dba, hsconsensus_CTCFfit)
newDBA_analysis<-dba.analyze(newDBA)

png("plots/025_diffbindMA_reimported_CTCFfit.png")
dba.plotMA(newDBA_analysis)
dev.off()
