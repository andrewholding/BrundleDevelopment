
### It will be like using counts, but we will use reads in peaks
library(DiffBind)

setwd("/Volumes/FlyPeaks/FlyPeaks")


# a is treated, b is not treated
#conversion<-setNames(
#  c("4b","2b","2a","1a","input","3b","3a","4a","1b"),
#  c("SLX-12998.D704_D505.C81G5ANXX.s_1.r_1.fq","SLX-12998.D704_D506.C81G5ANXX.s_1.r_1.fq","SLX-12998.D704_D507.C81G5ANXX.s_1.r_1.fq","SLX-12998.D705_D506.C81G5ANXX.s_1.r_1.fq","SLX-12998.D705_D507.C81G5ANXX.s_1.r_1.fq","SLX-12998.D705_D508.C81G5ANXX.s_1.r_1.fq","SLX-12998.D706_D505.C81G5ANXX.s_1.r_1.fq","SLX-12998.D706_D507.C81G5ANXX.s_1.r_1.fq","SLX-12998.D706_D508.C81G5ANXX.s_1.r_1.fq")
#)


# Read consensus peaks
hsconsensus<-readRDS("Rdata/015_SLX-12998_hsconsensus.rds")
mmconsensus<-readRDS("Rdata/015_SLX-12998_mmconsensus.rds")

# Get counts
hscounts<-hsconsensus[,-c(1:3)]
mmcounts<-mmconsensus[,-c(1:3)]
head(hscounts)
head(mmcounts)

### MA RPM in peaks
hsrpm<-apply(hscounts,2,function(x){
  1E6*x/sum(x)
})
mmrpm<-apply(mmcounts,2,function(x){
  1E6*x/sum(x)
})
M<-apply(hsrpm,1,function(x){
  fulvestrant<-mean(x[c(1,2,4,5)])
  untreated<-mean(x[c(3,6,7,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A<-apply(hsrpm,1,function(x){
  return(log10(sum(x)))
})
png("plots/016_SLX-12998_MA_RPMpeaks.png",w=1000,h=1000,p=30)
plot(A,M,pch=20,xlab="A, log10(RPM)",ylab="M, log2FC(fulvestrant)",main="RPM reads in peaks")
abline(h=0)
dev.off()


### MA RPM aligned reads
load("Rdata/012_SLX-12998_aligned.rda")
aligned<-aligned[grep("SLX-12998.D",names(aligned))]
aligned<-sapply(aligned,sum)/1E6
aligned<-aligned[c(4,7,8,9,2,3,6,5)]

hsrpm<-hscounts
for(i in 1:length(hsrpm)){
  hsrpm[i]<-hscounts[i]/aligned[i]
}
mmrpm<-mmcounts
for(i in 1:length(mmrpm)){
  mmrpm[i]<-mmcounts[i]/aligned[i]
}
M<-apply(hsrpm,1,function(x){
  fulvestrant<-mean(x[c(1,2,4,5)])
  untreated<-mean(x[c(3,6,7,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A<-apply(hsrpm,1,function(x){
  return(log10(sum(x)))
})
png("plots/016_SLX-12998_MA_RPMaligned.png",w=1000,h=1000,p=30)
plot(A,M,pch=20,xlab="A, log10(RPM)",ylab="M, log2FC(fulvestrant)",main="RPM aligned reads")
abline(h=0)
dev.off()


### MA RPM total reads
load("Rdata/012_SLX-12998_nrreads.rda")
nrreads<-nrreads[grep("SLX-12998.D",names(nrreads))]
nrreads<-nrreads/1E6
nrreads<-nrreads[c(4,7,8,9,2,3,6,5)]
hsrpm<-hscounts
for(i in 1:length(hsrpm)){
  hsrpm[i]<-hscounts[i]/nrreads[i]
}
mmrpm<-mmcounts
for(i in 1:length(mmrpm)){
  mmrpm[i]<-mmcounts[i]/nrreads[i]
}
M<-apply(hsrpm,1,function(x){
  fulvestrant<-mean(x[c(1,2,4,5)])
  untreated<-mean(x[c(3,6,7,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A<-apply(hsrpm,1,function(x){
  return(log10(sum(x)))
})
png("plots/016_SLX-12998_MA_RPMtotal.png",w=1000,h=1000,p=30)
plot(A,M,pch=20,xlab="A, log10(RPM)",ylab="M, log2FC(fulvestrant)",main="RPM total reads")
abline(h=0)
dev.off()



### MA Counts human
# a is treated, b is not treated
M<-apply(hscounts,1,function(x){
  fulvestrant<-mean(x[c(1,2,4,5)])
  untreated<-mean(x[c(3,6,7,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A<-apply(hscounts,1,function(x){
  return(log10(sum(x)))
})
png("plots/016_SLX-12998_MA_counts.png",w=1000,h=1000,p=30)
plot(A,M,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)",main="Raw counts in peaks")
abline(h=0)
dev.off()



### MA for Mouse + Human Counts
# a is treated, b is not treated
Mhs<-apply(hscounts,1,function(x){
  fulvestrant<-mean(x[c(1,2,4,5)])
  untreated<-mean(x[c(3,6,7,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Ahs<-apply(hscounts,1,function(x){
  return(log10(sum(x)))
})
Mmm<-apply(mmcounts,1,function(x){
  fulvestrant<-mean(x[c(1,2,4,5)])
  untreated<-mean(x[c(3,6,7,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Amm<-apply(mmcounts,1,function(x){
  return(log10(sum(x)))
})
png("plots/016_SLX-12998_MA_counts_HsMm.png",w=1000,h=1000,p=30)
plot(Ahs,Mhs,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)",main="Raw counts in peaks",ylim=c(-6.25,4))
points(Amm,Mmm,pch=20,col="darkolivegreen3")
abline(h=0)
legend("topright",legend=c("Mouse","Human"),pch=20,col=c("darkolivegreen3","black"))
lm1<-lm(Mmm~Amm)
abline(lm1$coef,col="red")
dev.off()


### 20th century normalization
# Residuals of the Drosophila fit
lm1<-lm(Mmm~Amm)
intercept<-lm1$coef[1]
angularcoeff<-lm1$coef[2]
MhsFit<-Mhs-(Ahs*angularcoeff)-intercept
MmmFit<-Mmm-(Amm*angularcoeff)-intercept

png("plots/016_SLX-12998-MA_counts_HsMm_Fit.png",w=1000,h=1000,p=30)
plot(Ahs,MhsFit,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)",main="Counts normalized by Mouse Distribution",ylim=c(-16,2))
points(Amm,MmmFit,pch=20,col="darkolivegreen3")
abline(h=0)
legend("topright",legend=c("Mouse","Human"),pch=20,col=c("darkolivegreen3","black"))
lm1<-lm(MmmFit~Amm)
abline(lm1$coef,col="red")
dev.off()





### MA for Mouse + Human Counts
# a is treated, b is not treated
Mhs<-apply(hscounts,1,function(x){
  fulvestrant<-mean(x[c(1,2,4,5)])
  untreated<-mean(x[c(3,6,7,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Ahs<-apply(hscounts,1,function(x){
  return(log10(sum(x)))
})
Mmm<-apply(mmcounts,1,function(x){
  fulvestrant<-mean(x[c(1,2,4,5)])
  untreated<-mean(x[c(3,6,7,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Amm<-apply(mmcounts,1,function(x){
  return(log10(sum(x)))
})
png("plots/016_SLX-12998_MA_counts_HsMm_twin_fit.png",w=1000,h=1000,p=30)
plot(Ahs,Mhs,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)",main="Raw counts in peaks",ylim=c(-6.25,4))
points(Amm,Mmm,pch=20,col="darkolivegreen3")
abline(h=0)
legend("topright",legend=c("Mouse","Human"),pch=20,col=c("darkolivegreen3","black"))
lm1<-lm(Mmm~Amm)
abline(lm1$coef,col="red")
lm1<-lm(Mhs~Ahs)
abline(lm1$coef,col="purple")
legend("bottomright",legend=c("Mouse","Human"),pch=20,col=c("red","purple"))

dev.off()

