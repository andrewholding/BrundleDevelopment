cd SLX-14229_mmhs/peaks

cat  SLX-14229.D701_D504.HJJL7BBXX.s_8.r_1.fq.gz_peaks.xls SLX-14229.D702_D503.HJJL7BBXX.s_8.r_1.fq.gz_peaks.xls > CTCF_merged.bed
grep -v \# CTCF_merged.bed  | grep -v start > CTCF_merged_filtered.bed
sort -k1,1 -k2,2n CTCF_merged_filtered.bed > CTCF_merged_sorted_filtered.bed
bedtools merge -i CTCF_merged_sorted_filtered.bed  >CTCF_union.bed

cd human

cat  SLX-14229.D701_D504.HJJL7BBXX.s_8.r_1.fq.gz_peaks.xls SLX-14229.D702_D503.HJJL7BBXX.s_8.r_1.fq.gz_peaks.xls > CTCF_merged.bed
grep -v \# CTCF_merged.bed  | grep -v start > CTCF_merged_filtered.bed
sort -k1,1 -k2,2n CTCF_merged_filtered.bed > CTCF_merged_sorted_filtered.bed
bedtools merge -i CTCF_merged_sorted_filtered.bed  >CTCF_union.bed

cd ..
cd mouse
cat  SLX-14229.D701_D504.HJJL7BBXX.s_8.r_1.fq.gz_peaks.xls SLX-14229.D702_D503.HJJL7BBXX.s_8.r_1.fq.gz_peaks.xls > CTCF_merged.bed
grep -v \# CTCF_merged.bed  | grep -v start > CTCF_merged_filtered.bed
sort -k1,1 -k2,2n CTCF_merged_filtered.bed > CTCF_merged_sorted_filtered.bed
bedtools merge -i CTCF_merged_sorted_filtered.bed  >CTCF_union.bed

