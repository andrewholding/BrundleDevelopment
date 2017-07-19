#!/bin/bash

echo Run with bash

hsMinReads=`cut -f2 -d"," SLX-8047_dmhs_titration/blacklist_filtered/human/reads.csv | sort -n | head -1`
dmMinReads=`cut -f2 -d"," SLX-8047_dmhs_titration/blacklist_filtered/drosophila/reads.csv | sort -n | head -1`

echo "Min Human:      $hsMinReads"
echo "Min Drosophila: $dmMinReads"


cd SLX-8047_dmhs_titration/blacklist_filtered/drosophila

mkdir 100


for f in *.r_1.fq.gz.bam
do
	echo -ne "$f, " 
	#libSize=`samtools flagstat $f |grep "mapped (" | grep -o "^\d*\ " | tr -d '\n'`
	#echo -ne "$libSize, "
	#echo "scale=3; $dmMinReads/$libSize"  | bc | awk '{printf "%.3f\n", $0}' |tr -d '\n' 
	#scalePercent=`echo "scale=2; $dmMinReads/$libSize"  | bc | awk '{printf "%.3f\n", $0}' |tr -d '\n'` 
	#if [ $scalePercent = "1.000" ]
	#then 
	#	samtools view -b $f > 100/$f
	#else
	#	samtools view -s $scalePercent -b $f > 100/$f
	#fi

	#gshuf is just shuf on Linux
	cat <(samtools view -H $f) <(samtools view -q 1  $f | gshuf -n $dmMinReads) | samtools view -bS - | samtools sort - -o 100/$f 
	samtools index 100/$f
	echo " done"
done

cd ../../..
