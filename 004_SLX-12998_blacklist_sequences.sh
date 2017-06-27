### Remove the blacklist from the alignments for SLX-12298
bl=../blacklists/hg19-blacklist.bed

cd ./SLX-12998
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




