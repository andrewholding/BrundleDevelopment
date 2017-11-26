cd SLX-15090/peaks/
mkdir CTCF

for f in *.broadPeak
do
	echo $f
	bedtools intersect -a $f -b CTCF_union.bed > CTCF/$f 
done

mkdir H4K12ac 
for f in *.broadPeak
do
	echo $f
	bedtools subtract -b CTCF_union.bed -a $f -A > H4K12ac/$f 
done


