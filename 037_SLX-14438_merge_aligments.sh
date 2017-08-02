mkdir ./SLX-14438_merged

cd ./SLX-14438_merged


for bam in $(find ../SLX-14438_mmhs/blacklist_filtered/human/ -type f);
do
	
	bam2=${bam//14438/14229}	
	bam3=${bam2/HKN27BBXX/HJJL7BBXX}
	samtools merge `basename $bam`  $bam  $bam3
done

mkdir sorted

for bam in *bam
do
	samtools sort $bam sorted/$bam
done

cd sorted

for bam in *bam
do
	samtools index $bam
done



