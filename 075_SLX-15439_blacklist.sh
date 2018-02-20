
bl=../blacklists/hg19-blacklist.bed


cd SLX-15439
mkdir ./blacklist_filtered
for f in *.bam
do
	echo $f
	bedtools intersect -v -abam $f -b $bl  > ./blacklist_filtered/$f
done

### Re-index

cd ./blacklist_filtered
for f in *.bam
do
samtools index $f
done

cd ..



