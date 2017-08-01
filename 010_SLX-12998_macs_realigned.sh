### MACS peak caller

### Run macs on the blacklisted data

mkdir ./SLX-12998_mmhs/peaks
cd ./SLX-12998_mmhs/peaks
control=../blacklist_filtered/SLX-12998.D708_D502.HH772BBXX.s_2.r_1.fq.gz.bam

for bam in ../blacklist_filtered/*bam
do
root=`basename $bam .bam`
macs2 callpeak -t $bam -c $control -f BAM -n $root -g hs
done




