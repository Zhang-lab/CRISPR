
#Only for EVT open regions (to TSC) that overlap to EVT KO vs WT close DARs
intersectBed -a EVT_open.bed -b EVT_KO_vs_WT_close.bed -wa -wb |cut -f 1-3 |sort -u > EVT_open_overlap_EVT_KOvsWT_close.bed

computeMatrix reference-point -R EVT_open_overlap_EVT_KOvsWT_close.bed -S EVT_KO.bw EVT_WT.bw  TSC_WT.bw --referencePoint center -a 5000 -b 5000 -bs 100 --missingDataAsZero -p 24 -o matrix_EVT_KOnWT_TSCWT_overEVTopen_overlap_EVT_KOvsWT_close
plotHeatmap -m matrix_EVT_KOnWT_TSCWT_overEVTopen_overlap_EVT_KOvsWT_close -o EVT_KOnWT_TSCWT_overEVTopen_EVT_KOvsWT_close.pdf --heatmapHeight 36 --perGroup --heatmapWidth 11 -min 0  --colorList '#ffffff,orange,#000000'

