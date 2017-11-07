cd SLX-14438_merged/peaks
bedtools subtract -A -a CTCF_union.bed -b ../../motifAnalysis/bed/output.bed > CTCF_EREremoved.bed
