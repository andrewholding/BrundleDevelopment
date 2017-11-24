### Run it on the blacklisted data
cd ./SLX-15091/peaks
control=../blacklist_filtered/SLX-15091.D702_D507.HMTVKBBXX.s_8.bwa.homo_sapiens.bam

for bam in ../blacklist_filtered/*bam
do
	root=`basename $bam .bam`
	macs2 callpeak -t $bam -c $control -f BAM -n $root -g hs
done


