library(ChIPpeakAnno)
library(GenomicRanges)


m1 = read.table("bed/ERoverlap/Carrol-GenesDev.bed", sep="\t")
m2 = read.table("bed/ERoverlap/ER_45minsGains_001_hg18.bed", sep="\t")
m3 = read.table("bed/ERoverlap/LinER.txt", sep="\t")
m4 = read.table("bed/ERoverlap/Welboren.bed", sep="\t")
m5 = read.table("bed/ERoverlap/NA_peaks.narrowPeak", sep="\t") #Called with macs2 against EtOH control
m5<-m5[c(1:3)]


m1.r = BED2RangedData(m1)
m2.r = BED2RangedData(m2)
m3.r = BED2RangedData(m3)
m4.r = BED2RangedData(m4)
m5.r = BED2RangedData(m5)
#mkd = makeVennDiagram(list(m1.r,m2.r,m3.r,m4.r,m5.r),
#                      NameOfPeaks=c("Ross-Innes 2010", "Guertin 2017", "Lin 2007", "Welboren 2009","Ceschin 2011"), totalTest=200000)

mkd = makeVennDiagram(list(m1.r,m2.r,m4.r,m5.r),
                      NameOfPeaks=c("Ross-Innes 2010", "Guertin 2017", "Welboren 2009","Ceschin 2011"))
