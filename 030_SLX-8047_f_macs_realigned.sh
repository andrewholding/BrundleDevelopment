### MACS peak caller
### RUn it on the blacklisted data

cd SLX-8047_dmhs_titration/blacklist_filtered/drosophila


control=SLX-8047.D705_D507.C81G5ANXX.s_1.r_1.fq.gz.bam
for percent in 1 2 5 10 15 20 100
#for percent in 100
do
	cd $percent
	mkdir peaks
	cd peaks
	for bam in ../*bam
	do
		root=`basename $bam .bam`
		echo $percent
		macs2 callpeak -t $bam -c ../$control -f BAM -n $root -g dm
	done
	cd ..
	cd ..
done
cd ..



