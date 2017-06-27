### Remove the blacklist from the alignments for SLX-8047
bl=../blacklists/hg19-blacklist.bed

cd ./SLX-8047
for f in *.bam
do
	echo $f
	bedtools intersect -v -abam $f -b $bl  > blacklist_filtered/$f
done

### Re-index
cd ./blacklist_filtered
for f in *.bam
do
samtools index $f
done

cd ..





