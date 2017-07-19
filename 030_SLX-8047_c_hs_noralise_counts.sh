#!/bin/bash

echo Run with Bash

hsMinReads=`cut -f2 -d"," SLX-8047_dmhs_titration/blacklist_filtered/human/reads.csv | sort -n | head -1`
dmMinReads=`cut -f2 -d"," SLX-8047_dmhs_titration/blacklist_filtered/drosophila/reads.csv | sort -n | head -1`

echo "Min Human:      $hsMinReads"
echo "Min Drosophila: $dmMinReads"


cd SLX-8047_dmhs_titration/blacklist_filtered/human

mkdir 100


for f in *.r_1.fq.gz.bam
do
	echo -ne "$f, " 
	#gshuf is just shuf on Linux
	cat <(samtools view -H $f) <(samtools view -q 1  $f | gshuf -n $hsMinReads) | samtools view -bS - | samtools sort - -o 100/$f 
	samtools index 100/$f
	echo " done"
done

#Ensure min read depth is met 
cp ../SLX-8047.D705_D506.C81G5ANXX.s_1.r_1.fq.gz.bam .
cp ../SLX-8047.D705_D506.C81G5ANXX.s_1.r_1.fq.gz.bam.bai .

cd ../../..
