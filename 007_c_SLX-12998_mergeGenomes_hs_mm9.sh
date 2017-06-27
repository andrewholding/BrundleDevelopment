### We should merge the genomes

#human=../Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
human=../Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
#droso=../Drosophila_melanogaster/NCBI/build5.41/Sequence/WholeGenomeFasta/genome.fa
#droso=../Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa
mouse=../Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa


mkdir genomes/mmhs
cd ./genomes/mmhs
sed "s/^>/>hs_/" $human > tmp
sed "s/^>/>mm_/" $mouse > tmp2
cat tmp tmp2 > mmhs.fa
rm tmp tmp2


### Merge Annotations
human=../Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-current/Genes/genes.gtf
#droso=../Drosophila_melanogaster/UCSC/dm3/Annotation/Archives/archive-current/Genes/genes.gtf
mouse=../Mus_musculus/UCSC/mm9/Annotation/Archives/archive-current/Genes/genes.gtf

sed "s/^/hs_/" $human > tmp
sed "s/^/mm_/" $mouse > tmp2

cat tmp tmp2 > mmhs.gtf
rm tmp tmp2

