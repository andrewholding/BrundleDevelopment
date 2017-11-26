cd SLX-15090/peaks
cat SLX-15090.D702_D504.HMTVKBBXX.s_7.bwa.homo_sapiens_peaks.narrowPeak SLX-15090.D702_D501.HMTVKBBXX.s_7.bwa.homo_sapiens_peaks.narrowPeak > CTCF_merged.bed
grep -v \# CTCF_merged.bed  | grep -v start > CTCF_merged_filtered.bed
sort -k1,1 -k2,2n CTCF_merged_filtered.bed > CTCF_merged_sorted_filtered.bed
bedtools merge -i CTCF_merged_sorted_filtered.bed  >CTCF_union.bed


