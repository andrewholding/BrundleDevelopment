### Run it on the blacklisted data
cd ./SLX-15090/peaks
control=../blacklist_filtered/SLX-15090.D702_D503.HMTVKBBXX.s_7.bwa.homo_sapiens.bam

for bam in ../blacklist_filtered/*bam
do
	root=`basename $bam .bam`
	#macs2 callpeak -t $bam -c $control -f BAM -n $root -g hs
	macs2 callpeak -t $bam -c $control --broad -f BAM -n $root -g hs
done


