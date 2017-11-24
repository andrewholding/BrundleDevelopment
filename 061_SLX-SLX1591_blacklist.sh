### Remove the blacklist from the alignments for merged alignement
bl=../blacklists/hg19-blacklist.bed 

cd ./SLX-15091

mkdir ./blackist_filtered
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





