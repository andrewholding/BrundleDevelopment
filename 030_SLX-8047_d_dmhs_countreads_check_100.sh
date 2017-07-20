#!/bin/bash
echo Must be run with bash.
echo
echo Drosophila

cd SLX-8047_dmhs_titration/blacklist_filtered/drosophila/100
for f in *.r_1.fq.gz.bam
do
	echo -ne "$f, " 
	samtools flagstat $f |grep "mapped (" | grep -o "^\d*\ " 
done

echo Human:

cd ../../human/100
for f in *.r_1.fq.gz.bam
do
	echo -ne "$f, " 
	samtools flagstat $f |grep "mapped (" | grep -o "^\d*\ " 
done

cd ../../..

