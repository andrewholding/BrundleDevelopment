cd SLX-14438_merged/peaks/
mkdir CTCF

for f in *.narrowPeak
do
	echo $f
	bedtools intersect -a $f -b CTCF_union.bed  > CTCF/$f 
done

mkdir ER
for f in *.narrowPeak
do
	echo $f
	bedtools subtract -b CTCF_union.bed -a $f -A > ER/$f 
done


