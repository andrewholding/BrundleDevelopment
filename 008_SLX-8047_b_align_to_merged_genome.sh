### Align, Sort and Build indeces
mkdir ./SLX-8047_dmhs
cd ./SLX-8047_dmhs
genome=../genomes/dmhs/dmhs
for fq in ../SLX-8047/*fq.gz
do
root=`basename $fq`
bowtie2 -p 32 -x $genome -U $fq > tmp.sam \
&& samtools view -Sbh tmp.sam > tmp.bam \
&& samtools sort tmp.bam ${root}  \
&& samtools index ${root}.bam
done
rm tmp.sam tmp.bam







