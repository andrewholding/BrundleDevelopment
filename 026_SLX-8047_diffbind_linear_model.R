setwd("/Volumes/FlyPeaks/FlyPeaks")
library(DiffBind)


filename<-"Rdata/026_SLX-8047_dba_human.rda"
if(!file.exists(filename)){
  dba<-dba(sampleSheet = "samplesheet/samplesheet_SLX8047_hs.csv")
  save(dba,file=filename)
} else {
  load(filename)
}


filename<-"Rdata/026_SLX-8047_dba.counts_human.rda"
if(!file.exists(filename)){
  dba <- dba.count(dba, summits=200)
  dba <- dba.count(dba, peaks=NULL, score=DBA_SCORE_READS)
  save(dba,file=filename)
} else {
  load(filename)
}

dba_analysis<-dba.analyze(dba)
png("plots/026_diffbindMA_import.png")
dba.plotMA(dba_analysis)
dev.off()


filename<-"Rdata/026_SLX-8047_dba_drosophila.rda"
if(!file.exists(filename)){
  dba_dm<-dba(sampleSheet = "samplesheet/samplesheet_SLX8047_dm.csv")
  save(dba_dm,file=filename)
} else {
  load(filename)
}


filename<-"Rdata/026_SLX-8047_dba.counts_drosophila.rda"
if(!file.exists(filename)){
  dba_dm <- dba.count(dba_dm, summits=200)
  dba_dm <- dba.count(dba_dm, peaks=NULL, score=DBA_SCORE_READS)
  save(dba_dm,file=filename)
} else {
  load(filename)
}


dba_analysis_dm<-dba.analyze(dba_dm)
png("plots/026_diffbindMA_import_dm.png")
dba.plotMA(dba_analysis_dm)
dev.off()



hsconsensus <- dba.peakset(dba_analysis, bRetrieve = T, DataType = DBA_DATA_FRAME)
dmconsensus <- dba.peakset(dba_analysis_dm, bRetrieve = T, DataType = DBA_DATA_FRAME)
#Correct names to samplesheet and test on intial data
names(hsconsensus)<-c("CHR","START","END","1a","1b","2a","2b","3a","3b","4a","4b")
names(dsconsensus)<-c("CHR","START","END","1a","1b","2a","2b","3a","3b","4a","4b")


#Now lets try a linear model.

load("Rdata/012_SLX-8047_aligned.rda")
#Remove lost reads bam
aligned<-aligned[grep("SLX-8047.D",names(aligned))]
aligned<-sapply(aligned,sum)/1E6
#remove inputs/controls to match sample order here
aligned<-aligned[c(4,9,3,2,7,6,8,1)]

hscounts<-hsconsensus[,-c(1:3)]
dmcounts<-dmconsensus[,-c(1:3)]

hsrpm<-hscounts
for(i in 1:length(hsrpm)){
  hsrpm[i]<-hscounts[i]/aligned[i]
}

dmrpm<-dmcounts
for(i in 1:length(dmrpm)){
  dmrpm[i]<-dmrpm[i]/aligned[i]
}



M_RPM<-apply(hsrpm,1,function(x){
  untreated<-mean(x[c(2,4,6,8)])
  fulvestrant<-mean(x[c(1,3,5,7)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A_RPM<-apply(hsrpm,1,function(x){
  return(log10(sum(x)))
})

M_dm_RPM<-apply(dmrpm,1,function(x){
  untreated<-mean(x[c(2,4,6,8)])
  fulvestrant<-mean(x[c(1,3,5,7)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A_dm_RPM<-apply(dmrpm,1,function(x){
  return(log10(sum(x)))
})

png("plots/026_MAplot_HsDm.png")
plot(A_RPM,M_RPM,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="RPM aligned reads")
points(A_dm_RPM,M_dm_RPM,pch=20,col="cornflowerblue")
lm1<-lm(M_dm_RPM~A_dm_RPM)
abline(lm1$coef,col="red4")
abline(h=0)
dev.off()

untreated<-rowMeans(dmrpm[c(1,3,5,7)])
fulvestrant<-rowMeans(dmrpm[c(2,4,6,8)])

png("plots/026_Dm_peak_count_comparision_log.png")
plot(log2(untreated),log2(fulvestrant), pch=20,
     xlab="Log2(counts) in peak after treatment" ,  ylab="Log2(counts) in peak before treatment" ,
     main="Comparision of Log2(counts) in peaks for Drosophila")

lm1<-lm(log2(untreated)~  log2(fulvestrant))
abline(lm1$coef,col="red3")
abline(h=0)

intercept<-lm1$coef[1]
angularcoeff<-lm1$coef[2] #is > 1 which implies less fuvestrant cells

fulvestrant_log2fit<-2**intercept*fulvestrant**angularcoeff
points(log2(fulvestrant_log2fit),log2(untreated),pch=20, col="royalblue3" )
abline(c(lm1$coef),col="red3")
abline(h=0)
lm1<-lm(log2(untreated) ~log2(fulvestrant_log2fit))
abline(c(lm1$coef),col="purple")
abline(h=0)
legend("topleft",legend=c("Raw", "Normalised"),pch=20,col=c("black","royalblue3"))
dev.off()


#Check fit by MA plot

log2fit_coeff<-angularcoeff
log2fit_offset<-intercept

M_RPM_log_corrected<-apply(hsrpm,1,function(x){
  untreated<-2**log2fit_offset*mean(x[c(2,4,6,8)])**log2fit_coeff
  fulvestrant<-mean(x[c(1,3,5,7)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A_RPM_log_corrected<-apply(hsrpm,1,function(x){
  untreated<-2**log2fit_offset*mean(x[c(2,4,6,8)])**log2fit_coeff
  fulvestrant<-mean(x[c(1,3,5,7)])
  return(log10(sum(fulvestrant+untreated)))
})

M_dm_RPM_log_corrected<-apply(dmrpm,1,function(x){
  untreated<-2**log2fit_offset*mean(x[c(2,4,6,8)])**(log2fit_coeff)
  fulvestrant<-mean(x[c(1,3,5,7)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A_dm_RPM_log_corrected<-apply(dmrpm,1,function(x){
  untreated<-2**log2fit_offset*mean(x[c(2,4,6,8)])**log2fit_coeff
  fulvestrant<-mean(x[c(1,3,5,7)])
  return(log10(sum(fulvestrant+untreated)))
})


png("plots/026_MAplot_dm_log2fit.png")
plot(A_RPM_log_corrected,M_RPM_log_corrected,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="RPM aligned reads - log2 CTCF Fit")
points(A_dm_RPM_log_corrected,M_dm_RPM_log_corrected,pch=20,col="cornflowerblue")
lm1<-lm(M_dm_RPM_log_corrected~A_dm_RPM_log_corrected)
abline(lm1$coef,col="red4")
abline(h=0)
dev.off()

#Apply the log2 dmfit
hsconsensus_dm_fit<-hsconsensus
dmconsensus_dm_fit<-dmconsensus
#change to RPM, fit and then convert back. 

hsconsensus_dm_fit[c(5,7,9,11)]<-
    round(2**(log2fit_offset)*(hsconsensus[c(5,7,9,11)]/aligned[c(2,4,6,8)])**(log2fit_coeff)*aligned[c(2,4,6,8)])


dmconsensus_dm_fit[c(5,7,9,11)]<-
  round(2**(log2fit_offset)*(dmconsensus[c(5,7,9,11)]/aligned[c(2,4,6,8)])**(log2fit_coeff )*aligned[c(2,4,6,8)])

lm1<-lm(log2(rowMeans(dmconsensus_dm_fit[c(4,6,8,10)]))~
        log2(rowMeans(dmconsensus_dm_fit[c(5,7,9,11)])))
lm1$coefficients

newDBA <- DiffBind:::pv.resetCounts(dba, hsconsensus_dm_fit)
newDBA_analysis<-dba.analyze(newDBA)

png("plots/026_diffbindMA_reimported_dm_fit.png")
dba.plotMA(newDBA_analysis)
dev.off()
