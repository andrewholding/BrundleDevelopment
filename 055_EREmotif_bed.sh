#Take the motif sites we got from a genome wide search and expan them +500bp both directions.

cd bed
wget https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes --no-check-certificate

bedtools slop -i ../motifAnalysis/bed/output.bed -b 500 -g hg19.chrom.sizes > CTCFblacklist.bed
