### Align, Sort and Build indeces
mkdir ./SLX-14229_mmhs
cd ./SLX-14229_mmhs
genome=../genomes/mmhs/mmhs
for fq in ../SLX-14229/*fq.gz
do
root=`basename $fq`
bowtie2 -p 32 -x $genome -U $fq > tmp.sam \
&& samtools view -Sbh tmp.sam > tmp.bam \
&& samtools sort tmp.bam ${root}  \
&& samtools index ${root}.bam
done
rm tmp.sam tmp.bam







