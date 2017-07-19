#!/bin/bash
echo Must be run with bash.
echo Coping SLX-8047_dmhs may take some time...

cp -r SLX-8047_dmhs SLX-8047_dmhs_titration

rm SLX-8047_dmhs_titration/blacklist_filtered/drosophila/reads.csv
rm SLX-8047_dmhs_titration/blacklist_filtered/human/reads.csv

cd SLX-8047_dmhs_titration/blacklist_filtered/drosophila
for f in *.r_1.fq.gz.bam
do
	echo $f
	echo -ne "$f, " >>reads.csv
	samtools flagstat $f |grep "mapped (" | grep -o "^\d*\ " >>reads.csv
	#samtools view $f   |wc -l >> reads
	#samtools view $f  |cut -f10 |sort |uniq |wc -l>>reads
done

cd ../human
for f in *.r_1.fq.gz.bam
do
	echo $f
	echo -ne "$f, " >>reads.csv
	samtools flagstat $f |grep "mapped (" | grep -o "^\d*\ " >>reads.csv
	#samtools view $f   |wc -l >> reads
	#samtools view $f  |cut -f10 |sort |uniq |wc -l>>reads
done

cd ../../..

