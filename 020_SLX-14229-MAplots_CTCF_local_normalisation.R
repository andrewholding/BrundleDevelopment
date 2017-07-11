

setwd("/Volumes/FlyPeaks/FlyPeaks")


# Read consensus peaks
hsconsensus_ER<-readRDS("Rdata/018_SLX-14229_hsconsensus_ER.rds")
hsconsensus_CTCF<-readRDS("Rdata/018_SLX-14229_hsconsensus_CTCF.rds")
mmconsensus_ER<-readRDS("Rdata/018_SLX-14229_mmconsensus_ER.rds")
mmconsensus_CTCF<-readRDS("Rdata/018_SLX-14229_mmconsensus_CTCF.rds")


filename<-"Rdata/020_SLX-14229_hsconsensus_localnorm_ER.rda"
if(!file.exists(filename)){
normalisation_window<-2000 #take all CTCF peaks this many base pairs up and down stream
hsconsensus_localnorm_ER<-hsconsensus_ER

#Normalise all ER peaks to the CTCF peaks 1kb either side.

for (n in rownames(hsconsensus_ER)){
min_location<-hsconsensus_ER[n,2]-(normalisation_window/2)
max_location<-hsconsensus_ER[n,3]-(normalisation_window/2)

control_peaks <-
  hsconsensus_CTCF[
  hsconsensus_CTCF[,1]==hsconsensus_ER[n,1] &
  hsconsensus_CTCF[,2]<=min_location &
  hsconsensus_CTCF[,3]<=max_location,]

normalisation_factors<-colSums(control_peaks[,-c(1:3)])/mean(colSums(control_peaks[,-c(1:3)]))

hsconsensus_localnorm_ER[n,-c(1:3)]<-(hsconsensus_ER[n,-c(1:3)]/normalisation_factors)
}
save(hsconsensus_localnorm_ER,file=filename)
} else {
  load(filename)
}

filename<-"Rdata/020_SLX-14229_hsconsensus_localnorm_CTCF.rda"
if(!file.exists(filename)){
#Normalise all CTCF peaks to the CTCF peaks 1kb either side for consistancy
hsconsensus_localnorm_CTCF<-hsconsensus_CTCF
for (n in rownames(hsconsensus_CTCF)){
  min_location<-hsconsensus_CTCF[n,2]-(normalisation_window/2)
  max_location<-hsconsensus_CTCF[n,3]-(normalisation_window/2)
  
  control_peaks <-
    hsconsensus_CTCF[
      hsconsensus_CTCF[,1]==hsconsensus_ER[n,1] &
        hsconsensus_CTCF[,2]<=min_location &
        hsconsensus_CTCF[,3]<=max_location,]
  normalisation_factors<-colSums(control_peaks[,-c(1:3)])/mean(colSums(control_peaks[,-c(1:3)]))

    
  hsconsensus_localnorm_CTCF[n,-c(1:3)]<-(hsconsensus_CTCF[n,-c(1:3)]/normalisation_factors)
}
save(hsconsensus_localnorm_CTCF,file=filename)
} else {
  load(filename)
}
# Get counts same as previous code
hscounts_ER<-hsconsensus_localnorm_ER[,-c(1:3)]
hscounts_CTCF<-hsconsensus_localnorm_CTCF[,-c(1:3)]
hscounts_before_ER<-hsconsensus_ER[,-c(1:3)]
hscounts_before_CTCF<-hsconsensus_CTCF[,-c(1:3)]
#mmcounts_ER<-mmconsensus_localnorm_ER[,-c(1:3)]
#mmcounts_CTCF<-mmconsensus_localnorm_CTCF[,-c(1:3)]
head(hscounts_ER)
head(hscounts_CTCF)


### MA for Mouse + Human Counts
# a is treated, b is not treated
Mhs_ER<-apply(hscounts_ER,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])
  untreated<-mean(x[c(1,3,6)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Ahs_ER<-apply(hscounts_ER,1,function(x){
  return(log10(sum(x)))
})

Mhs_before_ER<-apply(hscounts_before_ER,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])
  untreated<-mean(x[c(1,3,6)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Ahs_before_ER<-apply(hscounts_before_ER,1,function(x){
  return(log10(sum(x)))
})

Mhs_before_CTCF<-apply(hscounts_before_CTCF,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])
  untreated<-mean(x[c(1,3,6)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Ahs_before_CTCF<-apply(hscounts_before_CTCF,1,function(x){
  return(log10(sum(x)))
})

#Mmm_ER<-apply(mmcounts_ER,1,function(x){
#  fulvestrant<-mean(x[c(2,4,5)])
#  untreated<-mean(x[c(1,3,6)])
#  fc<-mean(fulvestrant)/mean(untreated)
#  log2fc<-log2(fc)
#  return(log2fc)
#})
#Amm_ER<-apply(mmcounts_ER,1,function(x){
#  return(log10(sum(x)))
#})
Mhs_CTCF<-apply(hscounts_CTCF,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])
  untreated<-mean(x[c(1,3,6)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Ahs_CTCF<-apply(hscounts_CTCF,1,function(x){
  return(log10(sum(x)))
})
#Mmm_CTCF<-apply(mmcounts_CTCF,1,function(x){
#  fulvestrant<-mean(x[c(2,4,5)])
#  untreated<-mean(x[c(1,3,6)])
#  fc<-mean(fulvestrant)/mean(untreated)
#  log2fc<-log2(fc)
#  return(log2fc)
#})
#Amm_CTCF<-apply(mmcounts_CTCF,1,function(x){
#  return(log10(sum(x)))
#})


library(scales)

png("plots/020_SLX-14229_MA_counts_HsMm_LocalNorm_simple.png",w=1000,h=1000,p=30)
plot(Ahs_ER,Mhs_ER,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="Local normalised counts in peaks",col="black")
points(Ahs_CTCF,Mhs_CTCF,pch=20,col="gray")

abline(h=0)
legend("topright",legend=c("Human ER",  "Human CTCF"),pch=20,col=c("black","gray"), cex=0.5)
lm1<-lm(Mhs_CTCF~Ahs_CTCF)
abline(lm1$coef,col="blue")
lm1<-lm(Mhs_ER~Ahs_ER)
abline(lm1$coef,col="purple")
legend("bottomright",legend=c("Human ER", "Human CTCF"),pch=20,col=c("purple","blue"), cex=0.5)

dev.off()
png("plots/020_SLX-14229_MA_counts_HsMm_LocalNorm.png",w=1000,h=1000,p=30)
plot(Ahs_before_ER,Mhs_before_ER,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="Local normalised counts in peaks",col="blue")
points(Ahs_ER,Mhs_ER,pch=20,col="black")
points(Ahs_before_CTCF,Mhs_before_CTCF,pch=20,col="lightblue")
points(Ahs_CTCF,Mhs_CTCF,pch=20,col="gray")

abline(h=0)
legend("topright",legend=c("Human ER",  "Human CTCF","Human ER Raw", "Human CTCF Raw"),pch=20,col=c("black","gray","blue","lightblue"), cex=0.5)
lm1<-lm(Mhs_CTCF~Ahs_CTCF)
abline(lm1$coef,col="blue")
lm1<-lm(Mhs_ER~Ahs_ER)
abline(lm1$coef,col="purple")
legend("bottomright",legend=c("Human ER", "Human CTCF"),pch=20,col=c("purple","blue"), cex=0.5)

dev.off()

# Residuals  fit
lm1<-lm(Mhs_CTCF~Ahs_CTCF)
intercept<-lm1$coef[1]
angularcoeff<-lm1$coef[2]
MhsFit_ER<-Mhs_ER-(Ahs_ER*angularcoeff)-intercept
MhsFit_CTCF<-Mhs_CTCF-(Ahs_CTCF*angularcoeff)-intercept


png("plots/020_SLX-14229_MA_counts_HsMm_LocalNorm_LM.png",w=1000,h=1000,p=30)
plot(Ahs_ER,MhsFit_ER,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="Corrected local normalised counts in peaks")
points(Ahs_CTCF,MhsFit_CTCF,pch=20,col="gray")
abline(h=0)
legend("topright",legend=c("Human ER",  "Human CTCF"),pch=20,col=c("black","gray"), cex=0.5)
lm1<-lm(MhsFit_CTCF~Ahs_CTCF)
abline(lm1$coef,col="blue")
lm1<-lm(MhsFit_ER~Ahs_ER)
abline(lm1$coef,col="purple")
legend("bottomright",legend=c("Human ER", "Human CTCF"),pch=20,col=c("purple","blue"), cex=0.5)

dev.off()


