
# The BAMs and BEDs should be split.

### BAMs, split in two
cd ./SLX-14229_mmhs/blacklist_filtered

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


### Peaks BEDs and XLSs, split in two
cd ../../peaks

echo `pwd`

mkdir human
mkdir mouse

for bed in *peaks.bed
do
grep -v "mm_chr" $bed | sed s/hs_chr/chr/g > human/$bed
grep -v "hs_chr" $bed | sed s/mm_chr/chr/g > mouse/$bed
done

for xls in *peaks.xls
do
grep -v "mm_chr" $xls | sed s/hs_chr/chr/g > human/$xls
grep -v "hs_chr" $xls | sed s/mm_chr/chr/g > mouse/$xls
done

for narrow in *narrowPeak
do
grep -v "mm_chr" $narrow | sed s/hs_chr/chr/g > human/$narrow
grep -v "hs_chr" $narrow | sed s/mm_chr/chr/g > mouse/$narrow
done
