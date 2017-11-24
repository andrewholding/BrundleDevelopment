cd SLX-15091/peaks

cat SLX-15091.D702_D508.HMTVKBBXX.s_8.bwa.homo_sapiens_peaks.xls SLX-15091.D702_D505.HMTVKBBXX.s_8.bwa.homo_sapiens_peaks.xls > CTCF_merged.bed
grep -v \# CTCF_merged.bed  | grep -v start > CTCF_merged_filtered.bed
sort -k1,1 -k2,2n CTCF_merged_filtered.bed > CTCF_merged_sorted_filtered.bed
bedtools merge -i CTCF_merged_sorted_filtered.bed  >CTCF_union.bed


