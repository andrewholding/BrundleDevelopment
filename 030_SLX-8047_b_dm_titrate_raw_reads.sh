#!/bin/bash


#Raw reads are needed as taking altering the human fly ratio breaks everything. 
#Therefore just down sample s2 reads
echo Run with bash
#Reads not including input
#MCF 	75454626
#S2 	64258555
# MGA put it lower

cd SLX-8047_dmhs_titration/blacklist_filtered/drosophila


for percent in 1 2 5 10 15 20 30 40 50
do
	mkdir A$percent
	for f in *.r_1.fq.gz.bam
	do
		scalePercent=`echo "scale=2; $percent/85"  | bc | tr -d '\n'` 
		echo $scalePercent
		echo "samtools view -s 0$scalePercent -b $f  | samtools sort - -o A$percent/$f"
		samtools view -s 0$scalePercent -b $f  | samtools sort - -o A$percent/$f
		samtools index A$percent/$f
		echo "$f, $percent,  done"
	done
done
