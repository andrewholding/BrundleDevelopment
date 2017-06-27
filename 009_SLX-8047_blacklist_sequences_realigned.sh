### Remove the blacklist from the alignments for merged alignement
bl=../blacklists/dmhs-blacklist.bed


cd ./SLX-8047_dmhs

mkdir ./blacklist_filtered
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





