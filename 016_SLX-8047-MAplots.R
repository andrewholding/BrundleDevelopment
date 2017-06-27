
### It will be like using counts, but we will use reads in peaks
library(DiffBind)

setwd("/Volumes/Estrocycle and FlyPeaks/FlyPeaks")


# a is treated, b is not treated
conversion<-setNames(
  c("4b","2b","2a","1a","input","3b","3a","4a","1b"),
  c("SLX-8047.D704_D505.C81G5ANXX.s_1.r_1.fq","SLX-8047.D704_D506.C81G5ANXX.s_1.r_1.fq","SLX-8047.D704_D507.C81G5ANXX.s_1.r_1.fq","SLX-8047.D705_D506.C81G5ANXX.s_1.r_1.fq","SLX-8047.D705_D507.C81G5ANXX.s_1.r_1.fq","SLX-8047.D705_D508.C81G5ANXX.s_1.r_1.fq","SLX-8047.D706_D505.C81G5ANXX.s_1.r_1.fq","SLX-8047.D706_D507.C81G5ANXX.s_1.r_1.fq","SLX-8047.D706_D508.C81G5ANXX.s_1.r_1.fq")
)


# Read consensus peaks
hsconsensus<-readRDS("Rdata/015_hsconsensus.rds")
dmconsensus<-readRDS("Rdata/015_dmconsensus.rds")

# Get counts
hscounts<-hsconsensus[,-c(1:3)]
dmcounts<-dmconsensus[,-c(1:3)]
head(hscounts)
head(dmcounts)

### MA RPM in peaks
hsrpm<-apply(hscounts,2,function(x){
  1E6*x/sum(x)
})
dmrpm<-apply(dmcounts,2,function(x){
  1E6*x/sum(x)
})
M<-apply(hsrpm,1,function(x){
  fulvestrant<-mean(x[c(1,3,5,7)])
  untreated<-mean(x[c(2,4,6,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A<-apply(hsrpm,1,function(x){
  return(log10(sum(x)))
})
png("plots/016_MA_RPMpeaks.png",w=1000,h=1000,p=30)
plot(A,M,pch=20,xlab="A, log10(RPM)",ylab="M, log2FC(fulvestrant)",main="RPM reads in peaks")
abline(h=0)
dev.off()


### MA RPM aligned reads
load("Rdata/012_aligned.rda")
aligned<-aligned[grep("SLX-8047.D",names(aligned))]
aligned<-sapply(aligned,sum)/1E6
names(aligned)<-conversion[names(aligned)]
hsrpm<-hscounts
for(i in 1:length(hsrpm)){
  hsrpm[i]<-hscounts[i]/aligned[i]
}
dmrpm<-dmcounts
for(i in 1:length(dmrpm)){
  dmrpm[i]<-dmcounts[i]/aligned[i]
}
M<-apply(hsrpm,1,function(x){
  fulvestrant<-mean(x[c(1,3,5,7)])
  untreated<-mean(x[c(2,4,6,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A<-apply(hsrpm,1,function(x){
  return(log10(sum(x)))
})
png("plots/016_MA_RPMaligned.png",w=1000,h=1000,p=30)
plot(A,M,pch=20,xlab="A, log10(RPM)",ylab="M, log2FC(fulvestrant)",main="RPM aligned reads")
abline(h=0)
dev.off()


### MA RPM total reads
load("Rdata/012_nrreads.rda")
nrreads<-nrreads[grep("SLX-8047.D",names(nrreads))]
nrreads<-nrreads/1E6
names(nrreads)<-conversion[names(nrreads)]
hsrpm<-hscounts
for(i in 1:length(hsrpm)){
  hsrpm[i]<-hscounts[i]/nrreads[i]
}
dmrpm<-dmcounts
for(i in 1:length(dmrpm)){
  dmrpm[i]<-dmcounts[i]/nrreads[i]
}
M<-apply(hsrpm,1,function(x){
  fulvestrant<-mean(x[c(1,3,5,7)])
  untreated<-mean(x[c(2,4,6,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A<-apply(hsrpm,1,function(x){
  return(log10(sum(x)))
})
png("plots/016_MA_RPMtotal.png",w=1000,h=1000,p=30)
plot(A,M,pch=20,xlab="A, log10(RPM)",ylab="M, log2FC(fulvestrant)",main="RPM total reads")
abline(h=0)
dev.off()



### MA Counts human
# a is treated, b is not treated
M<-apply(hscounts,1,function(x){
  fulvestrant<-mean(x[c(1,3,5,7)])
  untreated<-mean(x[c(2,4,6,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
A<-apply(hscounts,1,function(x){
  return(log10(sum(x)))
})
png("plots/016_MA_counts.png",w=1000,h=1000,p=30)
plot(A,M,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)",main="Raw counts in peaks")
abline(h=0)
dev.off()



### MA for Drosophila + Human Counts
# a is treated, b is not treated
Mhs<-apply(hscounts,1,function(x){
  fulvestrant<-mean(x[c(1,3,5,7)])
  untreated<-mean(x[c(2,4,6,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Ahs<-apply(hscounts,1,function(x){
  return(log10(sum(x)))
})
Mdm<-apply(dmcounts,1,function(x){
  fulvestrant<-mean(x[c(1,3,5,7)])
  untreated<-mean(x[c(2,4,6,8)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Adm<-apply(dmcounts,1,function(x){
  return(log10(sum(x)))
})
png("plots/016_MA_counts_HsDm.png",w=1000,h=1000,p=30)
plot(Ahs,Mhs,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)",main="Raw counts in peaks",ylim=c(-6.25,2))
points(Adm,Mdm,pch=20,col="cornflowerblue")
abline(h=0)
legend("topright",legend=c("Drosophila","Human"),pch=20,col=c("cornflowerblue","black"))
lm1<-lm(Mdm~Adm)
abline(lm1$coef,col="red")
dev.off()


### 20th century normalization
# Residuals of the Drosophila fit
lm1<-lm(Mdm~Adm)
intercept<-lm1$coef[1]
angularcoeff<-lm1$coef[2]
MhsFit<-Mhs-(Ahs*angularcoeff)-intercept
MdmFit<-Mdm-(Adm*angularcoeff)-intercept

png("plots/016_MA_counts_HsDm_Fit.png",w=1000,h=1000,p=30)
plot(Ahs,MhsFit,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)",main="Counts normalized by Drosophila Distribution",ylim=c(-6.25,2))
points(Adm,MdmFit,pch=20,col="cornflowerblue")
abline(h=0)
legend("topright",legend=c("Drosophila","Human"),pch=20,col=c("cornflowerblue","black"))
lm1<-lm(MdmFit~Adm)
abline(lm1$coef,col="red")
dev.off()






