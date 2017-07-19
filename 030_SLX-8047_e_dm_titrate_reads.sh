#!/bin/bash

echo Run with bash

hsMinReads=`cut -f2 -d"," SLX-8047_dmhs_titration/blacklist_filtered/human/reads.csv | sort -n | head -1`
dmMinReads=`cut -f2 -d"," SLX-8047_dmhs_titration/blacklist_filtered/drosophila/reads.csv | sort -n | head -1`

echo "Min Human:      $hsMinReads"
echo "Min Drosophila: $dmMinReads"


cd SLX-8047_dmhs_titration/blacklist_filtered/drosophila


for percent in 1 2 5 10 15 20
do
	mkdir $percent
	for f in *.r_1.fq.gz.bam
	do
		scalePercent=`echo "$hsMinReads*$percent/100"  | bc | tr -d '\n'` 
		cat <(samtools view -H $f) <(samtools view -q 1  $f | gshuf -n $scalePercent) | samtools view -bS - | samtools sort - -o $percent/$f
		samtools index $percent/$f
		echo "$f, $percent,  done"
	done
done
cd ../../..
