cd SLX-14229_mmhs/peaks/

cd human

cat CTCF_merged_sorted_filtered.bed | awk '$7>20' > CTCF_merged_sorted_filtered_cutoff.bed
bedtools merge -i CTCF_merged_sorted_filtered_cutoff.bed  >CTCF_union_cutoff.bed

mkdir CTCF2

for f in *.narrowPeak
do
	echo $f
	bedtools intersect -a $f -b CTCF_union_cutoff.bed  > CTCF2/$f 
done

mkdir ER2
for f in *.narrowPeak
do
	echo $f
	bedtools subtract -b CTCF_union_cutoff.bed -a $f -A > ER2/$f 
done
