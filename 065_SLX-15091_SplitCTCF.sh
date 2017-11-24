cd SLX-15091/peaks/
mkdir CTCF

for f in *.narrowPeak
do
	echo $f
	bedtools intersect -a $f -b CTCF_ERConsensusremoved.bed > CTCF/$f 
done

mkdir ER
for f in *.narrowPeak
do
	echo $f
	bedtools subtract -b CTCF_union.bed -a $f -A > ER/$f 
done


