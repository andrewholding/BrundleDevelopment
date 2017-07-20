### MACS peak caller
### RUn it on the blacklisted data

cd SLX-8047_dmhs_titration/blacklist_filtered/human/100


control=SLX-8047.D705_D507.C81G5ANXX.s_1.r_1.fq.gz.bam
mkdir peaks
cd peaks
pwd
for bam in ../*bam
do
	root=`basename $bam .bam`
	echo $percent
	macs2 callpeak -t $bam -c ../$control -f BAM -n $root -g hs
done




