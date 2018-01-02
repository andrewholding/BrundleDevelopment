#cat ~/homer/data/genomes/hg19/hg19.tss  |cut -f2-4  > tss.bed to generate TSS.bed
cd bed
bedtools intersect -wa -a ERConsensus500.bed -b tss.bed  > ERTSS.bed

cd ../SLX-15090/peaks/H4K12ac
mkdir EROnlyTSS
cd EROnlyTSS

#Only take peaks near ER using bed/ERTSS.bed we generated earlier

for file in ../*.broadPeak
do
	 bedtools intersect -wa -a $file -b ../../../../bed/ERTSS.bed > `basename $file`
done

