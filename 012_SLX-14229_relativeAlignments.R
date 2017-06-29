# TODO: Add comment
# 
# Author: giorgi01 & holdin01
###############################################################################

##### Load everything
## Number of reads per sample from data from CRUK pipeline before realignment.
setwd("/Volumes/Flypeaks/flypeaks/")
filename<-"Rdata/012_SLX-14229_nrreads.rda"
if(!file.exists(filename)){
  lines<-system("zgrep -c $ SLX-14229/*fq.gz",intern=TRUE)
  tmp<-strsplit(lines,":")
  nrreads<-sapply(tmp,function(x){as.numeric(x[2])/4})
  names(nrreads)<-sapply(tmp,function(x){x[1]})
  names(nrreads)<-gsub("SLX-14229/","",names(nrreads))
  names(nrreads)<-gsub(".HJJL7BBXX.s_8.r_1.fq.gz","",names(nrreads))
  names(nrreads)<-gsub(".r_1.fq.gz","",names(nrreads))
  save(nrreads,file=filename)
} else {
  load(filename)
}

## Reads aligned (on genomic portions too) on realigned data

filename<-"Rdata/012_SLX-14229_aligned.rda"

if(!file.exists(filename)){
  aligned<-list()
  for(nr in names(nrreads)){
    bn<-paste0(nr,".HJJL7BBXX.s_8.r_1.fq.gz.bam")
    message(bn)
    file<-paste0("./SLX-14229_mmhs/",bn)
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


nrreads<-nrreads[1:10]#Remove lost reads
aligned<-aligned[1:10]

### Prepare relative alignments
samples<-names(nrreads)
toplot<-matrix(nrow=length(samples),ncol=5)
colnames(toplot)<-c("MmTot","HsTot","Tot","MmRel","HsRel")
rownames(toplot)<-samples

mms<-sapply(aligned,function(x){
  sum(x[grep("mm_",names(x))])
})
hss<-sapply(aligned,function(x){
  sum(x[grep("hs_",names(x))])
})
toplot[,1]<-mms
toplot[,2]<-hss
toplot[,3]<-mms+hss
toplot[,4]<-mms/(mms+hss)
toplot[,5]<-hss/(mms+hss)

toplot<-toplot[grep('SLX-14229.D', sort(rownames(toplot)), value=TRUE),]

## Barplot of relative alignments
# Convert sample names to informative ones

rownames(toplot)<-c(
  "1+",
  "1-",
  "CTCF-",
  "Input-",
  "2+",
  "CTCF+",
  "2-",
  "3-",
  "Input+",
  "3+"
)



# Sort by total reads aligned
toplot<-toplot[order(-toplot[,3]),]

png("plots/012_SLX-14229_relativeAligned.png",w=1000,h=1000,point=30)
par(las=2,mar=c(4,4,3,5))
bp<-barplot(t(toplot[,4:5]),beside=TRUE,main="Relative Read Alignment in samples",
            ylab="Fraction of total aligned",col=c("royalblue2","tomato"),ylim=c(-0.05,1.05))
grid()
legend("topright",legend=c("Mouse","Human"),col=c("royalblue2","tomato"),pch=15)
bpx<-apply(bp,2,mean)

raxis<-(toplot[,3])/(max(toplot[,3]))
lines(bpx,raxis,type="b",lwd=3,pch=16)
axis(4,at=seq(0,1,0.2),
     labels=seq(0,1,0.2) * (max(toplot[,3]))
     )
par(las=3)
mtext("Reads aligned",side=4,line=4)

dev.off()























