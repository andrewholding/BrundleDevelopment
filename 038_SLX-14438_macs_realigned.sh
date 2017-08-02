### MACS peak caller

### RUn it on the blacklisted data
mkdir ./SLX-14438_merged/peaks
cd ./SLX-14438_merged/peaks
control=../sorted/SLX-14438.D703_D503.HKN27BBXX.s_8.r_1.fq.gz.bam.bam

for bam in ../sorted/*bam
do
root=`basename $bam .bam`
macs2 callpeak -t $bam -c $control -f BAM -n $root -g hs
done




