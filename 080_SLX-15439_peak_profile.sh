cd SLX-15439/blacklist_filtered/
makeTagDirectory PDX01/ SLX-15439.D701_D501.HNKNYBBXX.s_6.bwa.homo_sapiens.bam
makeTagDirectory PDX03/ SLX-15439.D701_D503.HNKNYBBXX.s_6.bwa.homo_sapiens.bam
makeTagDirectory PDX04/ SLX-15439.D702_D501.HNKNYBBXX.s_6.bwa.homo_sapiens.bam
makeTagDirectory PDX05/ SLX-15439.D702_D503.HNKNYBBXX.s_6.bwa.homo_sapiens.bam

annotatePeaks.pl ../../txt//CTCF_ERConsensusremoved.bed hg19 -size 6000 -hist 25 -d PDX01 PDX03 PDX04 PDX05 > Homer_CTCF_Profile.txt
annotatePeaks.pl ../../txt/PDXconsensus.txt hg19 -size 6000 -hist 25 -d PDX01 PDX03 PDX04 PDX05 > Homer_ER_Profile.txt
