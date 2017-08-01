


setwd("/Volumes/FlyPeaks/FlyPeaks")


# Read consensus peaks
hsconsensus_ER<-readRDS("Rdata/018_SLX-14229_hsconsensus_ER.rds")
hsconsensus_CTCF<-readRDS("Rdata/018_SLX-14229_hsconsensus_CTCF.rds")
mmconsensus_ER<-readRDS("Rdata/018_SLX-14229_mmconsensus_ER.rds")
mmconsensus_CTCF<-readRDS("Rdata/018_SLX-14229_mmconsensus_CTCF.rds")


filename<-"Rdata/020_SLX-14229_hsconsensus_localnorm_ER.rda"
load(filename)

filename<-"Rdata/020_SLX-14229_hsconsensus_localnorm_CTCF.rda"
load(filename)
}

# Get counts same as previous code
hscounts_localnorm_ER<-hsconsensus_localnorm_ER[,-c(1:3)]
hscounts_localnorm_CTCF<-hsconsensus_localnorm_CTCF[,-c(1:3)]
hscounts_ER<-hsconsensus_ER[,-c(1:3)]
hscounts_CTCF<-hsconsensus_CTCF[,-c(1:3)]
head(hscounts_ER)
head(hscounts_CTCF)




### Mean vs sd for with and without local normalisation

Mean_localnorm_ER_fulvestant<-apply(hscounts_localnorm_ER,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])
  return(fulvestrant)
})

Mean_localnorm_ER_untreated<-apply(hscounts_localnorm_ER,1,function(x){
  untreated<-mean(x[c(1,3,6)])
  return(untreated)
})

Mean_ER_fulvestant<-apply(hscounts_ER,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])
  return(fulvestrant)
})

Mean_ER_untreated<-apply(hscounts_ER,1,function(x){
  untreated<-mean(x[c(1,3,6)])
  return(untreated)
})

#sd
sd_localnorm_ER_fulvestant<-apply(hscounts_localnorm_ER,1,function(x){
  fulvestrant<-sd(x[c(2,4,5)])
  return(fulvestrant)
})

sd_localnorm_ER_untreated<-apply(hscounts_localnorm_ER,1,function(x){
  untreated<-sd(x[c(1,3,6)])
  return(untreated)
})

sd_ER_fulvestant<-apply(hscounts_ER,1,function(x){
  fulvestrant<-sd(x[c(2,4,5)])
  return(fulvestrant)
})

sd_ER_untreated<-apply(hscounts_ER,1,function(x){
  untreated<-sd(x[c(1,3,6)])
  return(untreated)
})


#Do I want std_dev of the mean?

png("plots/021_SLX-14229_sd_vs_Mean.png",w=2000,h=1000,p=30)

par(mfrow=c(1,2))
plot(Mean_ER_fulvestant,sd_ER_fulvestant,pch=20,col="gray",xlab="Mean",ylab="SD", main="SD of treated peak counts")
points(Mean_localnorm_ER_fulvestant,sd_localnorm_ER_fulvestant,pch=20,col="goldenrod3")
lm1<-lm(sd_ER_fulvestant~Mean_ER_fulvestant)
abline(lm1$coef,col="grey")
lm1<-lm(sd_localnorm_ER_fulvestant~Mean_localnorm_ER_fulvestant)
abline(lm1$coef,col="goldenrod3")
legend("topleft",legend=c("Raw", "Normalised"),pch=20,col=c("grey","goldenrod3"), cex=0.5)


plot(Mean_ER_untreated,sd_ER_untreated,pch=20,col="gray",xlab="Mean",ylab="SD", main="SD of untreated peak counts")
points(Mean_localnorm_ER_untreated,sd_localnorm_ER_untreated,pch=20,col="goldenrod3")
lm1<-lm(sd_ER_untreated~Mean_ER_untreated)
abline(lm1$coef,col="grey")
lm1<-lm(sd_localnorm_ER_untreated~Mean_localnorm_ER_untreated)
abline(lm1$coef,col="goldenrod3")
legend("topleft",legend=c("Raw", "Normalised"),pch=20,col=c("grey","goldenrod3"), cex=0.5)

dev.off()

### MA RPM aligned reads
load("Rdata/012_SLX-12998_aligned.rda")
aligned<-aligned[grep("SLX-12998.D",names(aligned))]
aligned<-sapply(aligned,sum)/40E6
aligned<-aligned[c(4,7,8,9,2,3,6,5)]

hsrpm<-hscounts_ER
for(i in 1:length(hsrpm)){
  hsrpm[i]<-hscounts_ER[i]/aligned[i]
}

Mean_ER_rpm_fulvestant<-apply(hsrpm,1,function(x){
  fulvestrant<-mean(x[c(2,4,5)])
  return(fulvestrant)
})

Mean_ER_rpm_untreated<-apply(hsrpm,1,function(x){
  untreated<-mean(x[c(1,3,6)])
  return(untreated)
})


sd_ER__rpm_fulvestant<-apply(hsrpm,1,function(x){
  fulvestrant<-sd(x[c(2,4,5)])
  return(fulvestrant)
})

sd_ER__rpm_untreated<-apply(hsrpm,1,function(x){
  untreated<-sd(x[c(1,3,6)])
  return(untreated)
})

#Reads per RP40M
points(Mean_ER_rpm_fulvestant,sd_ER__rpm_fulvestant,pch=20,col="blue")
lm1<-lm(sd_ER_fulvestant~Mean_ER_fulvestant)
abline(lm1$coef,col="blue")

points(Mean_ER_rpm_untreatedt,sd_ER__rpm_untreated,pch=20,col="blue")
lm1<-lm(sd_ER_untreated~Mean_ER_untreated)
abline(lm1$coef,col="blue")

