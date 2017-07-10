cd SLX-14229_mmhs/peaks/mouse
mkdir CTCF

for f in *.narrowPeak
do
	echo $f
	bedtools intersect -a CTCF_union.bed -b $f > CTCF/$f 
done

mkdir ER
for f in *.narrowPeak
do
	echo $f
	bedtools subtract -b CTCF_union.bed -a $f -A > ER/$f 
done


cd ..
cd human

mkdir CTCF

for f in *.narrowPeak
do
	echo $f
	bedtools intersect -a CTCF_union.bed -b $f > CTCF/$f 
done

mkdir ER
for f in *.narrowPeak
do
	echo $f
	bedtools subtract -b CTCF_union.bed -a $f -A > ER/$f 
done
