# TODO: Add comment
# 
# Author: giorgi01
###############################################################################

##ANH: My bowtie isn't patched for GZIP so had to unzip first.

### Align, Sort and Build indeces
mkdir ./SLX-12998_mmhs
cd ./SLX-12998_mmhs
genome=../genomes/mmhs/mmhs
for fq in ../SLX-12998/*fq.gz
do
root=`basename $fq`
bowtie2 -p 32 -x $genome -U $fq > tmp.sam \
&& samtools view -Sbh tmp.sam > tmp.bam \
&& samtools sort tmp.bam ${root}  \
&& samtools index ${root}.bam
done
rm tmp.sam tmp.bam







