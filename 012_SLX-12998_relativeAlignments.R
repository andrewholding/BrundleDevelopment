# TODO: Add comment
# 
# Author: giorgi01 & holdin01
###############################################################################

##### Load everything
## Number of reads per sample from data from CRUK pipeline before realignment.
setwd("/Volumes/Flypeaks/flypeaks/")
filename<-"Rdata/012_SLX-12998_nrreads.rda"
if(!file.exists(filename)){
  lines<-system("zgrep -c $ SLX-12998/*fq.gz",intern=TRUE)
  tmp<-strsplit(lines,":")
  nrreads<-sapply(tmp,function(x){as.numeric(x[2])/4})
  names(nrreads)<-sapply(tmp,function(x){x[1]})
  names(nrreads)<-gsub("SLX-12998/","",names(nrreads))
  names(nrreads)<-gsub(".H772BBXX.s_2.r_1.bwa.homo_sapiens.fq.gz","",names(nrreads))
  names(nrreads)<-gsub(".r_1.fq.gz","",names(nrreads))
  save(nrreads,file=filename)
} else {
  load(filename)
}

## Reads aligned (on genomic portions too) on realigned data

filename<-"Rdata/012_SLX-12998_aligned.rda"

if(!file.exists(filename)){
  aligned<-list()
  for(nr in names(nrreads)){
    bn<-paste0(nr,".r_1.fq.gz.bam")
    file<-paste0("./SLX-12998_mmhs/",bn)
    # SAM tools flags
    # -F256 removes the non primary alignments
    # -F4 remove the non aligned ones
    lines<-system(paste0("samtools view -F 260 ",file," | cut -f3 | sort | uniq -c"),intern=TRUE)
    tmp<-strsplit(lines," ")
    al<-as.numeric(sapply(tmp,function(x){
      return(x[length(x)-1])
    }))
    names(al)<-sapply(tmp,function(x){
      return(x[length(x)])
    })
    aligned[[nr]]<-al
  }
  save(aligned,file=filename)
} else {
  load(filename)
}

#NEEDS CHANGE FOR SLX-12998
nrreads<-nrreads[2:10]#Remove lost reads
aligned<-aligned[2:10]

### Prepare relative alignments
samples<-names(nrreads)
toplot<-matrix(nrow=length(samples),ncol=5)
colnames(toplot)<-c("DmTot","HsTot","Tot","DmRel","HsRel")
rownames(toplot)<-samples

dms<-sapply(aligned,function(x){
  sum(x[grep("dm_",names(x))])
})
hss<-sapply(aligned,function(x){
  sum(x[grep("hs_",names(x))])
})
toplot[,1]<-dms
toplot[,2]<-hss
toplot[,3]<-dms+hss
toplot[,4]<-dms/(dms+hss)
toplot[,5]<-hss/(dms+hss)

toplot<-toplot[grep('SLX-12998.D', sort(rownames(toplot)), value=TRUE),]

## Barplot of relative alignments
# Convert sample names to informative ones

rownames(toplot)<-c(
  "4b-",#"SLX-12998.D704_D505",
  "2b-",#"SLX-12998.D704_D506",
  "2a+",#"SLX-12998.D704_D507",
  "1a+",#"SLX-.D705_D506",
  "Input",#"SLX-.D705_D507",
  "3b-",#"SLX-.D705_D508",
  "3a+",#"SLX-12998.D706_D505",
  "4a+",#"SLX-12998.D706_D507",
  "1b-"#"SLX-12998.D706_D508"
)



# Sort by total reads aligned
toplot<-toplot[order(-toplot[,3]),]

png("plots/012_SLX-12998_relativeAligned.png",w=1000,h=1000,point=30)
par(las=2,mar=c(4,4,3,5))
bp<-barplot(t(toplot[,4:5]),beside=TRUE,main="Relative Read Alignment in samples",
            ylab="Fraction of total aligned",col=c("darkolivegreen4","tomato"),ylim=c(-0.05,1.05))
grid()
legend("topright",legend=c("Drosophila","Human"),col=c("darkolivegreen4","tomato"),pch=15)
bpx<-apply(bp,2,mean)

raxis<-(toplot[,3]-min(toplot[,3]))/(max(toplot[,3])-min(toplot[,3]))
lines(bpx,raxis,type="b",lwd=3,pch=16)
axis(4,at=pretty(raxis),
     labels=(pretty(raxis) * (max(toplot[,3])-min(toplot[,3])) + min(toplot[,3])
     ))
par(las=3)
mtext("Reads aligned",side=4,line=4)

dev.off()























