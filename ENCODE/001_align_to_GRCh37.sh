# TODO: Add comment
# 
# Author: A N Holding
###############################################################################

for bam in *bam
do
bamToFastq -i $bam -fq $bam.fq
done








