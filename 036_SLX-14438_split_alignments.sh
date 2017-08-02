
# The BAMs and BEDs should be splitted in Drosophila (for normalization) and in HUman (for diffbind)



### BAMs, split in two
cd ./SLX-14438_mmhs/blacklist_filtered

mkdir human
mkdir mouse

for bam in *bam
do
samtools view -h $bam | grep -v "mm_chr" | sed s/hs_chr/chr/g | samtools view -bS - > human/$bam
samtools view -h $bam | grep -v "hs_chr" | sed s/mm_chr/chr/g | samtools view -bS - > mouse/$bam 
done

# Index
cd ./human
echo `pwd`

for bam in *bam
do
samtools index $bam
done


cd ../mouse
echo `pwd`

for bam in *bam
do
samtools index $bam
done


