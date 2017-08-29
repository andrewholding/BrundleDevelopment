cd csv
awk -F":" '$1=$1' OFS="\t" 029_SLX14229_controlResultsDeseq.csv | tr -s ',' '\t' | tr -d '"' | awk -F"\t" -v OFS="\t" '{sub(/-/,"\t",$2)}1' | awk '$9 != "NA"' |tail -n +2 |less |awk '$9 < 0.01' > 044_significantCTCFPeaks.bed
cd ..
mkdir motifAnalysis
cd motifAnalysis
findMotifsGenome.pl ../csv/044_significantCTCFPeaks.bed hg19 ctcfBinding

mkdir bed
cd bed
scanMotifGenomeWide.pl ../ctcfBinding/homerMotifs.all.motifs hg19 -bed -int -keepAll > output.bed

