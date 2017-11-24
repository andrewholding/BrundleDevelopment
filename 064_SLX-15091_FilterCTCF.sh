#cd ./bed"
#bedtools slop -i ERConsensus.bed -b 500 -g hg19.chrom.sizes > ERConsensus500.bed
cd ./SLX-15091/peaks
bedtools subtract -A -a CTCF_union.bed -b ../../bed/ERConsensus500.bed > CTCF_ERConsensusremoved.bed
