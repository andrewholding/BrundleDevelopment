cd SLX-14438_merged/peaks

cat  SLX-14438.D701_D504.HKN27BBXX.s_8.r_1.fq.gz.bam_peaks.xls SLX-14438.D702_D503.HKN27BBXX.s_8.r_1.fq.gz.bam_peaks.xls > CTCF_merged.bed
grep -v \# CTCF_merged.bed  | grep -v start > CTCF_merged_filtered.bed
sort -k1,1 -k2,2n CTCF_merged_filtered.bed > CTCF_merged_sorted_filtered.bed
bedtools merge -i CTCF_merged_sorted_filtered.bed  >CTCF_union.bed


