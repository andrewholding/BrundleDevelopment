# TODO: Add comment
# 
# Author: A N Holding
###############################################################################


genome=../genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
for fq in *fq
do
bowtie2 -p 3 -x $genome -U $fq > tmp.sam \
&& samtools view -bh tmp.sam > tmp.bam \
&& samtools sort -o $fq.hg19.bam tmp.bam \
&& samtools index $fq.hg19.bam
done
rm tmp.sam tmp.bam







