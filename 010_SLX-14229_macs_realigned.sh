### MACS peak caller

### Run macs on the blacklisted data

mkdir ./SLX-14229_mmhs/peaks
cd ./SLX-14229_mmhs/peaks
control=../blacklist_filtered/SLX-14229.D703_D503.HJJL7BBXX.s_8.r_1.fq.gz.bam

for bam in ../blacklist_filtered/*bam
do
root=`basename $bam .bam`
macs2 callpeak -t $bam -c $control -f BAM -n $root -g hs
done




