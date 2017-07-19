### MACS peak caller
### RUn it on the blacklisted data

mkdir peaks
cd peaks
control=SLX-8047.D705_D507.C81G5ANXX.s_1.r_1.fq.gz.bam
for $percent in 1 2 5 10 15 20 
	for bam in ../*bam
	do
		root=`basename $bam .bam`
		echo "macs2 callpeak -t $bam -c ../$control -f BAM -n $percent -g dm"
	done
done
cd ..



