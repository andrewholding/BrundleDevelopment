cd SLX-15090/peaks/H4K12ac
mkdir EROnly
cd EROnly

#Only take peaks near ER using bed/ERConsensus.bed we generated earlier

for file in ../*.broadPeak
do
	 bedtools intersect -wa -a $file -b ../../../../bed/ERConsensus.bed > `basename $file`
done

