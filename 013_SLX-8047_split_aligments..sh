
# The BAMs and BEDs should be splitted in Drosophila (for normalization) and in HUman (for diffbind)



### BAMs, split in two
cd ./SLX-8047_dmhs/blacklist_filtered

mkdir human
mkdir drosophila

for bam in *bam
do
samtools view -h $bam | grep -v "dm_chr" | sed s/hs_chr/chr/g | samtools view -bS - > human/$bam
samtools view -h $bam | grep -v "hs_chr" | sed s/dm_chr/chr/g | samtools view -bS - > drosophila/$bam 
done

# Index
cd ./human
echo `pwd`

for bam in *bam
do
samtools index $bam
done


cd ../drosophila
echo `pwd`

for bam in *bam
do
samtools index $bam
done


### Peaks BEDs and XLSs, split in two
cd ../../peaks

echo `pwd`

mkdir human
mkdir drosophila

for bed in *peaks.bed
do
grep -v "dm_chr" $bed | sed s/hs_chr/chr/g > human/$bed
grep -v "hs_chr" $bed | sed s/dm_chr/chr/g > drosophila/$bed
done

for xls in *peaks.xls
do
grep -v "dm_chr" $xls | sed s/hs_chr/chr/g > human/$xls
grep -v "hs_chr" $xls | sed s/dm_chr/chr/g > drosophila/$xls
done

for narrow in *narrowPeak
do
grep -v "dm_chr" $narrow | sed s/hs_chr/chr/g > human/$narrow
grep -v "hs_chr" $narrow | sed s/dm_chr/chr/g > drosophila/$narrow
done
