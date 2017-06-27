### MACS peak caller

### RUn it on the blacklisted data
cd ./SLX-14229/peaks

#Used the I without ICI added
control=../blacklist_filtered/SLX-14229.D703_D503.HJJL7BBXX.s_8.bwa.homo_sapiens.bam

for bam in ../blacklist_filtered/*bam
do
root=`basename $bam .bam`
macs2 callpeak -t $bam -c $control -f BAM -n $root -g hs
done




