
### It will be like using counts, but we will use reads in peaks
library(DiffBind)

setwd("/Volumes/FlyPeaks/FlyPeaks")


# a is treated, b is not treated
conversion<-setNames(
  c("X4b","X2b","X2a","X1a","input","X3b","X3a","X4a","X1b"),
  c("SLX-8047.D704_D505.C81G5ANXX.s_1","SLX-8047.D704_D506.C81G5ANXX.s_1","SLX-8047.D704_D507.C81G5ANXX.s_1","SLX-8047.D705_D506.C81G5ANXX.s_1","SLX-8047.D705_D507.C81G5ANXX.s_1","SLX-8047.D705_D508.C81G5ANXX.s_1","SLX-8047.D706_D505.C81G5ANXX.s_1","SLX-8047.D706_D507.C81G5ANXX.s_1","SLX-8047.D706_D508.C81G5ANXX.s_1")
)


# Read consensus peaks
hsconsensus<-readRDS("Rdata/015_SLX-8047_hsconsensus.rds")
dmconsensus<-readRDS("Rdata/015_SLX-8047_dmconsensus.rds")

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
png("plots/016_SLX-8047_MA_RPMpeaks.png")
plot(A,M,pch=20,ylab=expression("log"[2]~"ChIP fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),main="RPM reads in peaks")
abline(h=0)
dev.off()


### MA RPM aligned reads
load("Rdata/012_SLX-8047_aligned.rda")
aligned<-aligned[grep("SLX-8047.D",names(aligned))]
aligned<-sapply(aligned,sum)/1E6
names(aligned)<-conversion[names(aligned)]
aligned = aligned[c(4,9,3,2,6,7,8,1)]
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
png("plots/016_SLX-8047_MA_RPMaligned.png")
plot(A,M,pch=20,ylab=expression("log"[2]~"ChIP fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),main="RPM aligned reads")
abline(h=0)
dev.off()


### MA RPM total reads
load("Rdata/012_SLX-8047_nrreads.rda")
nrreads<-nrreads[grep("SLX-8047.D",names(nrreads))]
nrreads<-nrreads/1E6
names(nrreads)<-conversion[names(nrreads)]
nrreads = nrreads[c(4,9,3,2,6,7,8,1)]
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
png("plots/016_SLX-8047_MA_RPMtotal.png")
plot(A,M,pch=20,ylab=expression("log"[2]~"ChIP fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),main="RPM total reads")
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
png("plots/016_SLX-8047_MA_counts.png")
plot(A,M,pch=20,ylab=expression("log"[2]~"ChIP fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),main="Raw counts in peaks")
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
png("plots/016_SLX-8047_MA_counts_HsDm.png")
plot(Ahs,Mhs,pch=20,ylab=expression("log"[2]~"ChIP fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),main="Raw counts in peaks",ylim=c(-6.25,2))
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

png("plots/016_SLX-8047-MA_counts_HsDm_Fit.png")
plot(Ahs,MhsFit,pch=20,ylab=expression("log"[2]~"ChIP fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),main="Counts normalized by Drosophila Distribution",ylim=c(-6.25,2))
points(Adm,MdmFit,pch=20,col="cornflowerblue")
abline(h=0)
legend("topright",legend=c("Drosophila","Human"),pch=20,col=c("cornflowerblue","black"))
lm1<-lm(MdmFit~Adm)
abline(lm1$coef,col="red")
dev.off()






