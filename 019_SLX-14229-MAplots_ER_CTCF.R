### It will be like using counts, but we will use reads in peaks
library(DiffBind)

setwd("/Volumes/FlyPeaks/FlyPeaks")


# a is treated, b is not treated
#conversion<-setNames(
#  c("4b","2b","2a","1a","input","3b","3a","4a","1b"),
#  c("SLX-14229.D704_D505.C81G5ANXX.s_1.r_1.fq","SLX-14229.D704_D506.C81G5ANXX.s_1.r_1.fq","SLX-14229.D704_D507.C81G5ANXX.s_1.r_1.fq","SLX-14229.D705_D506.C81G5ANXX.s_1.r_1.fq","SLX-14229.D705_D507.C81G5ANXX.s_1.r_1.fq","SLX-14229.D705_D508.C81G5ANXX.s_1.r_1.fq","SLX-14229.D706_D505.C81G5ANXX.s_1.r_1.fq","SLX-14229.D706_D507.C81G5ANXX.s_1.r_1.fq","SLX-14229.D706_D508.C81G5ANXX.s_1.r_1.fq")
#)


# Read consensus peaks
hsconsensus_ER<-readRDS("Rdata/018_SLX-14229_hsconsensus_ER.rds")
hsconsensus_CTCF<-readRDS("Rdata/018_SLX-14229_hsconsensus_CTCF.rds")
mmconsensus_ER<-readRDS("Rdata/018_SLX-14229_mmconsensus_ER.rds")
mmconsensus_CTCF<-readRDS("Rdata/018_SLX-14229_mmconsensus_CTCF.rds")


# Get counts
hscounts_ER<-hsconsensus_ER[,-c(1:3)]
hscounts_CTCF<-hsconsensus_CTCF[,-c(1:3)]
mmcounts_ER<-mmconsensus_ER[,-c(1:3)]
mmcounts_CTCF<-mmconsensus_CTCF[,-c(1:3)]
head(hscounts_ER)
head(hscounts_CTCF)


### MA for Mouse + Human Counts
# a is treated, b is not treated
Mhs_ER<-apply(hscounts_ER,1,function(x){
  fulvestrant<-mean(x[c(2,4,6)])
  untreated<-mean(x[c(1,3,5)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Ahs_ER<-apply(hscounts_ER,1,function(x){
  return(log10(sum(x)))
})
Mmm_ER<-apply(mmcounts_ER,1,function(x){
  fulvestrant<-mean(x[c(2,4,6)])
  untreated<-mean(x[c(1,3,5)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Amm_ER<-apply(mmcounts_ER,1,function(x){
  return(log10(sum(x)))
})
Mhs_CTCF<-apply(hscounts_CTCF,1,function(x){
  fulvestrant<-mean(x[c(2,4,6)])
  untreated<-mean(x[c(1,3,5)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Ahs_CTCF<-apply(hscounts_CTCF,1,function(x){
  return(log10(sum(x)))
})
Mmm_CTCF<-apply(mmcounts_CTCF,1,function(x){
  fulvestrant<-mean(x[c(2,4,6)])
  untreated<-mean(x[c(1,3,5)])
  fc<-mean(fulvestrant)/mean(untreated)
  log2fc<-log2(fc)
  return(log2fc)
})
Amm_CTCF<-apply(mmcounts_CTCF,1,function(x){
  return(log10(sum(x)))
})


library(scales)

png("plots/019_SLX-14229_MA_counts_HsMm.png",w=1000,h=1000,p=30)
plot(Ahs_ER,Mhs_ER,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="Raw counts in peaks")
points(Ahs_CTCF,Mhs_CTCF,pch=20,col="gray")
points(Amm_ER,Mmm_ER,pch=20,col="darkolivegreen3")
points(Amm_CTCF,Mmm_CTCF,pch=20,col="darkolivegreen")
abline(h=0)
legend("topright",legend=c("Mouse ER","Human ER", "Mouse CTCF", "Human CTCF"),pch=20,col=c("darkolivegreen3","black","darkolivegreen3","gray"), cex=0.5)
lm1<-lm(Mmm_CTCF~Amm_CTCF)
abline(lm1$coef,col="red")
lm1<-lm(Mhs_CTCF~Ahs_CTCF)
abline(lm1$coef,col="blue")
lm1<-lm(Mmm_ER~Amm_ER)
abline(lm1$coef,col="orange")
lm1<-lm(Mhs_ER~Ahs_ER)
abline(lm1$coef,col="purple")
legend("bottomright",legend=c("Mouse ER","Human ER", "Mouse CTCF", "Human CTCF"),pch=20,col=c("orange","purple","red","blue"), cex=0.5)

dev.off()


### 20th century normalization
# Residuals  fit
lm1<-lm(Mhs_CTCF~Ahs_CTCF)
intercept<-lm1$coef[1]
angularcoeff<-lm1$coef[2]
MhsFit_ER<-Mhs_ER-(Ahs_ER*angularcoeff)-intercept
MmmFit_ER<-Mmm_ER-(Amm_ER*angularcoeff)-intercept
MhsFit_CTCF<-Mhs_CTCF-(Ahs_CTCF*angularcoeff)-intercept
MmmFit_CTCF<-Mmm_CTCF-(Amm_CTCF*angularcoeff)-intercept



png("plots/019_SLX-14229_MA_counts_HsMm_Normalised.png",w=1000,h=1000,p=30)
plot(Ahs_ER,MhsFit_ER,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="Raw counts in peaks")
points(Ahs_CTCF,MhsFit_CTCF,pch=20,col="gray")
points(Amm_ER,MmmFit_ER,pch=20,col="darkolivegreen3")
points(Amm_CTCF,MmmFit_CTCF,pch=20,col="darkolivegreen")
abline(h=0)
legend("topright",legend=c("Mouse ER","Human ER", "Mouse CTCF", "Human CTCF"),pch=20,col=c("darkolivegreen3","black","darkolivegreen3","gray"), cex=0.5)
lm1<-lm(MmmFit_CTCF~Amm_CTCF)
abline(lm1$coef,col="red")
lm1<-lm(MhsFit_CTCF~Ahs_CTCF)
abline(lm1$coef,col="blue")
lm1<-lm(MmmFit_ER~Amm_ER)
abline(lm1$coef,col="orange")
lm1<-lm(MhsFit_ER~Ahs_ER)
abline(lm1$coef,col="purple")
legend("bottomright",legend=c("Mouse ER","Human ER", "Mouse CTCF", "Human CTCF"),pch=20,col=c("orange","purple","red","blue"), cex=0.5)

dev.off()

png("plots/019_SLX-14229_MA_counts_HsMm_Normalised_Alpha.png",w=1000,h=1000,p=30)

plot(Ahs_ER,MhsFit_ER,pch=20,xlab="A, log10(counts)",ylab="M, log2FC(fulvestrant)", main="Normalised to Human CTCF",col=alpha("black",0.2))
points(Ahs_CTCF,MhsFit_CTCF,pch=20,col=alpha("gray",0.2))
points(Amm_ER,MmmFit_ER,pch=20,col=alpha("darkolivegreen3",0.2))
points(Amm_CTCF,MmmFit_CTCF,pch=20,col=alpha("darkolivegreen",0.2))
abline(h=0)
legend("topright",legend=c("Mouse ER","Human ER", "Mouse CTCF", "Human CTCF"),pch=20,col=c("darkolivegreen3","black","darkolivegreen3","gray"), cex=0.5)
lm1<-lm(MmmFit_CTCF~Amm_CTCF)
abline(lm1$coef,col="red")
lm1<-lm(MhsFit_CTCF~Ahs_CTCF)
abline(lm1$coef,col="blue")
lm1<-lm(MmmFit_ER~Amm_ER)
abline(lm1$coef,col="orange")
lm1<-lm(MhsFit_ER~Ahs_ER)
abline(lm1$coef,col="purple")
legend("bottomright",legend=c("Mouse ER","Human ER", "Mouse CTCF", "Human CTCF"),pch=20,col=c("orange","purple","red","blue"), cex=0.5)
dev.off()




