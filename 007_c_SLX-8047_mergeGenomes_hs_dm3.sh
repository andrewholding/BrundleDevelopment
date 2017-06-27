### We should merge the genomes

#human=../Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
human=../Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
#droso=../Drosophila_melanogaster/NCBI/build5.41/Sequence/WholeGenomeFasta/genome.fa
droso=../Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa


mkdir genomes/dmhs
cd ./genomes/dmhs
sed "s/^>/>hs_/" $human > tmp
sed "s/^>/>dm_/" $droso > tmp2
cat tmp tmp2 > dmhs.fa
rm tmp tmp2


### Merge Annotations
human=../Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-current/Genes/genes.gtf
droso=../Drosophila_melanogaster/UCSC/dm3/Annotation/Archives/archive-current/Genes/genes.gtf

sed "s/^/hs_/" $human > tmp
sed "s/^/dm_/" $droso > tmp2

cat tmp tmp2 > dmhs.gtf
rm tmp tmp2

