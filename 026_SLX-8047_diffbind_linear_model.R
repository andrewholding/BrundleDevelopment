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



#use a stringent overlap for control peaks.
filename<-"Rdata/026_SLX-8047_dba.counts_drosophila.rda"
if(!file.exists(filename)){
  dba_dm <- dba.count(dba_dm, summits=200, minOverlap = 5)
  dba_dm <- dba.count(dba_dm, peaks=NULL, score=DBA_SCORE_READS)
  save(dba_dm,file=filename)
} else {
  load(filename)
}


dba_analysis_dm<-dba.analyze(dba_dm)
png("plots/026_diffbindMA_import_dm.png")
dba.plotMA(dba_analysis_dm)
dev.off()



hsconsensus <- dba.peakset(dba, bRetrieve = T, DataType = DBA_DATA_FRAME)
dmconsensus <- dba.peakset(dba_dm, bRetrieve = T, DataType = DBA_DATA_FRAME)
#Correct names to samplesheet and test on intial data
names(hsconsensus)<-c("CHR","START","END","1a","1b","2a","2b","3a","3b","4a","4b")
names(dmconsensus)<-c("CHR","START","END","1a","1b","2a","2b","3a","3b","4a","4b")


#Now lets try to correct the bias model.

load("Rdata/012_SLX-8047_aligned.rda")
#Remove lost reads bam
aligned<-aligned[grep("SLX-8047.D",names(aligned))]
aligned<-sapply(aligned,sum)/1E6
#remove inputs/controls to match sample order here
aligned<-aligned[c(4,9,3,2,7,6,8,1)]

#Corrected for Hs/Dm only Bams. Future work is to add code to do automatically
hs_aligned<-c(4222714,6996508,11021281,10222291,10098809,5988566,4303772,10222291)
hs_aligned<-sapply(hs_aligned,sum)/1E6
dm_aligned<-c(12457121,10702574,9040006,8234849,6539216,5527308,16722010,6231542 )
dm_aligned<-sapply(dm_aligned,sum)/1E6

aligned=hs_aligned+dm_aligned


hscounts<-hsconsensus[,-c(1:3)]
dmcounts<-dmconsensus[,-c(1:3)]



M<-apply(hscounts,1,function(x){
  untreated<-mean(x[c(2,4,6,8)])
  fulvestrant<-mean(x[c(1,3,5,7)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A<-apply(hscounts,1,function(x){
  return(log10(sum(x)))
})

M_dm<-apply(dmcounts,1,function(x){
  untreated<-mean(x[c(2,4,6,8)])
  fulvestrant<-mean(x[c(1,3,5,7)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A_dm<-apply(dmcounts,1,function(x){
  return(log10(sum(x)))
})

png("plots/026_MAplot_HsDm.png")
plot(A,M,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="Raw aligned reads")
points(A_dm,M_dm,pch=20,col="cornflowerblue")
lm1<-lm(M_dm~A_dm)
abline(lm1$coef,col="red4")
abline(h=0)
dev.off()

untreated<-rowMeans(dmcounts[c(1,3,5,7)])
fulvestrant<-rowMeans(dmcounts[c(2,4,6,8)])

png("plots/026_Dm_peak_count_comparision.png")
plot(untreated,fulvestrant, pch=20,
     xlab="Counts in peak after treatment" ,  ylab="Counts in peak before treatment" ,
     main="Comparision of Counts in peaks for Drosophila")
abline(0,1,col="grey")
lm1<-lm(  fulvestrant~0+ untreated)
abline(c(0,lm1$coef),col="red3")

intercept<-0
angularcoeff<-1/lm1$coef[1] #is < 1 which implies less fuvestrant cells

points(fulvestrant*angularcoeff+intercept,untreated,pch=20, col="royalblue3" )
fulvestrant_fit<-fulvestrant*angularcoeff+intercept
lm1<-lm(fulvestrant_fit  ~0+untreated)
abline(c(0,lm1$coef),col="purple")
legend("topleft",legend=c("Raw", "Normalised"),pch=20,col=c("black","royalblue3"))
dev.off()


#Check fit by MA plot

fit_coeff<-angularcoeff
fit_offset<-intercept

M_corrected<-apply(hscounts,1,function(x){
  untreated<-fit_coeff*mean(x[c(2,4,6,8)])+fit_offset
  fulvestrant<-mean(x[c(1,3,5,7)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A_corrected<-apply(hscounts,1,function(x){
  untreated<-fit_coeff*mean(x[c(2,4,6,8)])+fit_offset
  fulvestrant<-mean(x[c(1,3,5,7)])
  return(log10(sum(fulvestrant+untreated)))
})

M_dm_corrected<-apply(dmcounts,1,function(x){
  untreated<-fit_coeff*mean(x[c(2,4,6,8)])+fit_offset
  fulvestrant<-mean(x[c(1,3,5,7)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A_dm_corrected<-apply(dmcounts,1,function(x){
  untreated<-fit_coeff*mean(x[c(2,4,6,8)])+fit_offset
  fulvestrant<-mean(x[c(1,3,5,7)])
  return(log10(sum(fulvestrant+untreated)))
})


png("plots/026_MAplot_dm_fit.png")
plot(A_corrected,M_corrected,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="Raw aligned reads - log2 Dm Fit")
points(A_dm_corrected,M_dm_corrected,pch=20,col="cornflowerblue")
lm1<-lm(M_dm_corrected~A_dm_corrected)
abline(lm1$coef,col="red4")
abline(h=0)
dev.off()

#Apply the  dmfit
hsconsensus_dm_fit<-hsconsensus
dmconsensus_dm_fit<-dmconsensus


#DeSeq called by Diffbind will divide though by library size, 
#therefore undo our work. Solution is to correct for that before
#hand with this control factor.
control_factor<-sum(aligned[c(2,4,6,8)])/sum(aligned[c(1,3,5,7)])

hsconsensus_dm_fit[c(5,7,9,11)]<-
    (fit_coeff*hsconsensus[c(5,7,9,11)]+fit_offset)/control_factor


dmconsensus_dm_fit[c(5,7,9,11)]<-
  (fit_coeff*dmconsensus[c(5,7,9,11)]+fit_offset)



newDBA <- DiffBind:::pv.resetCounts(dba, hsconsensus_dm_fit)
newDBA_analysis<-dba.analyze(newDBA)


png("plots/026_diffbindMA_reimported_dm_fit.png")
dba.plotMA(newDBA_analysis)
dev.off()

