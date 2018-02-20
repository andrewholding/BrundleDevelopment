cd SLX-15439/peaks/
mkdir CTCF

#BEDs for CTCF from SLX-15091

for f in *.narrowPeak
do
	echo $f
	bedtools intersect -a $f -b ../../txt/CTCF_ERConsensusremoved.bed > CTCF/$f 
done

mkdir ER
for f in *.narrowPeak
do
	echo $f
	bedtools subtract -b ../../txt/CTCF_union.bed -a $f -A > ER/$f 
done

