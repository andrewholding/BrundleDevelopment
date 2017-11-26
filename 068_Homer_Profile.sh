cd SLX-15090/
makeTagDirectory H4K12ac-E2-ChIP-Seq/ SLX-15090.D702_D502.HMTVKBBXX.s_7.bwa.homo_sapiens.bam SLX-15090.D703_D501.HMTVKBBXX.s_7.bwa.homo_sapiens.bam SLX-15090.D702_D502.HMTVKBBXX.s_7.bwa.homo_sapiens.bam 
makeTagDirectory H4K12ac-None-ChIP-Seq/  SLX-15090.D701_D501.HMTVKBBXX.s_7.bwa.homo_sapiens.bam  SLX-15090.D701_D504.HMTVKBBXX.s_7.bwa.homo_sapiens.bam  SLX-15090.D703_D502.HMTVKBBXX.s_7.bwa.homo_sapiens.bam

annotatePeaks.pl tss hg19 -size 6000 -hist 25 -d H4K12ac-E2-ChIP-Seq/ H4K12ac-None-ChIP-Seq/ > Homer_Profile.txt
annotatePeaks.pl ../bed/ERConsensus.bed hg19 -size 6000 -hist 25 -d H4K12ac-E2-ChIP-Seq/ H4K12ac-None-ChIP-Seq/ > Homer_ER_Profile.txt
annotatePeaks.pl ../motifAnalysis/bed/output.bed hg19 -size 6000 -hist 25 -d H4K12ac-E2-ChIP-Seq/ H4K12ac-None-ChIP-Seq/ > Homer_ERE_Profile.txt


