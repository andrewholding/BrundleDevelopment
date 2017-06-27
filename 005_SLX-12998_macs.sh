### MACS peak caller
# It uses python!
# Downloaded version 1.4.2 from http://liulab.dfci.harvard.edu/MACS/Download.html
# This doesn't work anymore: sudo dpkg -i macs_1.4.2.deb

# So...
#cd /Users/giorgi01/Programs
#wget https://github.com/downloads/taoliu/MACS/MACS-1.4.2-1.tar.gz
#python setup.py install --user

# We need to install also PeakSplitter 1.0 (which I did from here: http://www.ebi.ac.uk/research/bertone/software)

#S2
#ICI
#1a +  D506 D705
#1b -  D508 D706
#2a +  D507 D704
#2b -  D506 D704
#3a +  D505 D706
#3b -  D508 D705
#4a +  D507 D706
#4b -  D505 D704
#Input D507 D705

#HC11
# "D707-D504" "ICI_4"
# "D708-D502" "Input Control"
# "D707-D507" "ICI_1"
# "D707-D508" "Control_4"
# "D708-D501" "ICI_3"
# "D707-D502" "ICI_2"
# "D707-D503" "Control_1"
# "D707-D501" "Input ICI"
# "D707-D506" "Control_3"
# "D707-D505" "Control_2"

# It runs now
# -t is the treatment bam file
# -c should be the control (Input D507 D705)
# -n is the experiment name
# -g is the genome size (Drosophila + Human = 3.2e+9, the size of dmhs.fa, but we will use the default "hs" for human genome)
# -w saves extended fragment pileup
# --call-subpeaks will use tghe PeakSplitter algorithm (requires -w)

### RUn it on the blacklisted data
cd ./SLX-12998/peaks
control=../blacklist_filtered/SLX-12998.D708_D502.HH772BBXX.s_2.bwa.homo_sapiens.bam

for bam in ../blacklist_filtered/*bam
do
root=`basename $bam .bam`
macs2 callpeak -t $bam -c $control -f BAM -n $root -g hs
done




