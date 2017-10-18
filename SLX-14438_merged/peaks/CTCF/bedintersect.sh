for bed in *.narrowPeak
do 
	bedtools intersect -a $bed -b ../../../bed/CTCF_union_removedER.bed > safe/$bed 
done

