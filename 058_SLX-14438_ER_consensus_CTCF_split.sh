cd SLX-14438_merged/peaks/
mkdir CTCF_noER

for f in *.narrowPeak
do
	echo $f
	bedtools intersect -a $f -b CTCF_ERConsensusremoved.bed  > CTCF_noER/$f 
done

