### MACS peak caller

### Run macs on the blacklisted data
cd ./SLX-15439
mkdir ./peaks
cd peaks


control=D702_D504
bam=D701_D501
macs2 callpeak -t ../blacklist_filtered/SLX-15439.$bam.HNKNYBBXX.s_6.bwa.homo_sapiens.bam -c ../blacklist_filtered/SLX-15439.$control.HNKNYBBXX.s_6.bwa.homo_sapiens.bam -f BAM -n $bam.bam -g hs &


control=D701_D504
bam=D703_D502
macs2 callpeak -t ../blacklist_filtered/SLX-15439.$bam.HNKNYBBXX.s_6.bwa.homo_sapiens.bam -c ../blacklist_filtered/SLX-15439.$control.HNKNYBBXX.s_6.bwa.homo_sapiens.bam -f BAM -n $bam.bam -g hs &


control=D702_D502
bam=D701_D503
macs2 callpeak -t ../blacklist_filtered/SLX-15439.$bam.HNKNYBBXX.s_6.bwa.homo_sapiens.bam -c ../blacklist_filtered/SLX-15439.$control.HNKNYBBXX.s_6.bwa.homo_sapiens.bam -f BAM -n $bam.bam -g hs &

control=D701_D502
bam=D702_D501
macs2 callpeak -t ../blacklist_filtered/SLX-15439.$bam.HNKNYBBXX.s_6.bwa.homo_sapiens.bam -c ../blacklist_filtered/SLX-15439.$control.HNKNYBBXX.s_6.bwa.homo_sapiens.bam -f BAM -n $bam.bam -g hs &

control=D703_D501
bam=D702_D503
macs2 callpeak -t ../blacklist_filtered/SLX-15439.$bam.HNKNYBBXX.s_6.bwa.homo_sapiens.bam -c ../blacklist_filtered/SLX-15439.$control.HNKNYBBXX.s_6.bwa.homo_sapiens.bam -f BAM -n $bam.bam -g hs &







