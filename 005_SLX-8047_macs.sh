### MACS peak caller
# -t is the treatment bam file
# -c should be the control (Input D507 D705)
# -n is the experiment name
# -g is the genome size (Drosophila + Human = 3.2e+9, the size of dmhs.fa, but we will use the default "hs" for human genome)
# -w saves extended fragment pileup
# --call-subpeaks will use tghe PeakSplitter algorithm (requires -w)

### Run it on the blacklisted data
cd ./SLX-8047/peaks
control=../blacklist_filtered/SLX-8047.D705_D507.C81G5ANXX.s_1.bwa.homo_sapiens.bam

for bam in ../blacklist_filtered/*bam
do
	root=`basename $bam .bam`
	macs2 callpeak -t $bam -c $control -f BAM -n $root -g hs
done




