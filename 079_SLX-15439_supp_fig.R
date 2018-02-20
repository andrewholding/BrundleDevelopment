library(Brundle)


pdx<-readRDS(file="Rdata/pdx.rds")
ctcf<-readRDS(file="Rdata/ctcf.rds")

pdx<- dba.count(pdx,peaks = NULL, score = DBA_SCORE_READS)
ctcf<- dba.count(ctcf,peaks = NULL, score = DBA_SCORE_READS)

peakset.ctcf <- dba.peakset(ctcf, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
peakset.er <- dba.peakset(pdx, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)

normalisation.factor<-colSums(peakset.ctcf[4:8])/median(colSums(peakset.ctcf[4:8]))



normalised.er<-cbind(peakset.er[1:3],peakset.er[4]/normalisation.factor[1],
                                     peakset.er[5]/normalisation.factor[2],
                                     peakset.er[6]/normalisation.factor[3],
                                     peakset.er[7]/normalisation.factor[4],
                                     peakset.er[8]/normalisation.factor[5])

colSums(normalised.er[4:8])
colSums(peakset.er[4:8])

normalised.ctcf<-cbind(peakset.ctcf[1:3],
                       peakset.ctcf[4]/normalisation.factor[1],
                       peakset.ctcf[5]/normalisation.factor[2],
                       peakset.ctcf[6]/normalisation.factor[3],
                       peakset.ctcf[7]/normalisation.factor[4],
                       peakset.ctcf[8]/normalisation.factor[5])

colSums(normalised.ctcf[4:8])

rpm<-as.numeric(pdx$class["Reads",])/10^6


addToDf<-function(data,factor,method,locus) {
        data.frame(counts=as.numeric(t(data)),factor=rep(factor,4),method=rep(method,4),rep=c("PDX01","PDX03","PDX04","PDX05"),locus=rep(locus,4))
}

library(ggpubr)

relativeRPM<-(rpm[c(1,3:5)]/median(rpm[c(1,3:5)]))

df<-addToDf(peakset.ctcf[13416,c(4,6:8)],"CTCF","Raw","RARa")
df<-rbind(df,addToDf(peakset.er[639,c(4,6:8)],"ER","Raw","RARa"))
df<-rbind(df,addToDf(peakset.ctcf[13416,c(4,6:8)]/relativeRPM,"CTCF","RPM","RARa"))
df<-rbind(df,addToDf(peakset.er[639,c(4,6:8)]/relativeRPM,"ER","RPM","RARa"))
df<-rbind(df,addToDf(normalised.ctcf[13416,c(4,6:8)],"CTCF","pfChIP","RARa"))
df<-rbind(df,addToDf(normalised.er[639,c(4,6:8)],"ER","pfChIP","RARa"))


ctcfSite<-"21382"
erSite<-"891"
name<-"GREB1"
df<-rbind(df,addToDf(peakset.ctcf[ctcfSite,c(4,6:8)],"CTCF","Raw",name))
df<-rbind(df,addToDf(peakset.er[erSite,c(4,6:8)],"ER","Raw",name))
df<-rbind(df,addToDf(normalised.ctcf[ctcfSite,c(4,6:8)],"CTCF","pfChIP",name))
df<-rbind(df,addToDf(normalised.er[erSite,c(4,6:8)],"ER","pfChIP",name))
df<-rbind(df,addToDf(peakset.ctcf[ctcfSite,c(4,6:8)]/relativeRPM,"CTCF","RPM",name))
df<-rbind(df,addToDf(peakset.er[erSite,c(4,6:8)]/relativeRPM,"ER","RPM",name))

ctcfSite<-"20453"
erSite<-"1064"
name<-"CLIC6"
df<-rbind(df,addToDf(peakset.ctcf[ctcfSite,c(4,6:8)],"CTCF","Raw",name))
df<-rbind(df,addToDf(peakset.er[erSite,c(4,6:8)],"ER","Raw",name))
df<-rbind(df,addToDf(normalised.ctcf[ctcfSite,c(4,6:8)],"CTCF","pfChIP",name))
df<-rbind(df,addToDf(normalised.er[erSite,c(4,6:8)],"ER","pfChIP",name))
df<-rbind(df,addToDf(peakset.ctcf[ctcfSite,c(4,6:8)]/relativeRPM,"CTCF","RPM",name))
df<-rbind(df,addToDf(peakset.er[erSite,c(4,6:8)]/relativeRPM,"ER","RPM",name))


ggbarplot(df, x="rep", y = "counts",
          fill = "method",    
          color = "white",            # Set bar border colors to white
          palette =c("#FC4E07", "#E7B800", "#00AFBB"),
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Corrected Counts",
          xlab = "PDX",
          legend.title = "Comparison",
          facet.by=c("locus","factor"),
          position=position_dodge(width=0.7)
)


ggboxplot(df, x="method", y = "counts",
          palette =c("#FC4E07", "#E7B800", "#00AFBB"),
          ylab = "Corrected Counts",
          xlab = "PDX",
          color="method",
          facet.by=c("locus","factor"),
          legend.title = "Comparison"
)
