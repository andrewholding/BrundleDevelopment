### blacklist DAC (for Human and Drosophila)s
# Empirically detected by ENCODE
# Downloaded from https://sites.google.com/site/anshulkundaje/projects/blacklists on August 30 2016
mkdir ./blacklists
cd ./blacklists
# Human hg19/GRCh37
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz --no-check-certificate  
gunzip wgEncodeDacMapabilityConsensusExcludable.bed.gz
mv wgEncodeDacMapabilityConsensusExcludable.bed hg19-blacklist.bed
# Mouse mm9
wget http://www.broadinstitute.org/~anshul/projects/mouse/blacklist/mm9-blacklist.bed.gz --no-check-certificate  
gunzip mm9-blacklist.bed.gz
# Nematode ce10
wget http://www.broadinstitute.org/~anshul/projects/worm/blacklist/ce10-blacklist.bed.gz --no-check-certificate  
gunzip ce10-blacklist.bed.gz
# Drosophila dm3
wget http://www.broadinstitute.org/~anshul/projects/fly/blacklist/dm3-blacklist.bed.gz --no-check-certificate  
gunzip dm3-blacklist.bed.gz


### Generate a merged Drosophila/Human blacklist
droso=./dm3-blacklist.bed
homo=./hg19-blacklist.bed

sed "s/^/hs_/" $homo |cut -f1-3 > tmp
sed "s/^/dm_/" $droso > tmp2
cat tmp tmp2 > dmhs-blacklist.bed
rm tmp tmp2

### Generate a merged Mouse/Human blacklist
mouse=./mm9-blacklist.bed

sed "s/^/hs_/" $homo |cut -f1-3 > tmp
sed "s/^/mm_/" $mouse > tmp2
cat tmp tmp2 > mmhs-blacklist.bed
rm tmp tmp2














