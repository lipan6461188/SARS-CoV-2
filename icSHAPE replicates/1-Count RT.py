
cd /Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing

icSHAPE-pipe countRT \
    -in 7.sam2tab/D1_virus.tab,7.sam2tab/D2_virus.tab,7.sam2tab/N1_virus.tab,7.sam2tab/N2_virus.tab,7.sam2tab/T1_virus.tab,7.sam2tab/T2_virus.tab \
    -size index/SARS2_complement.size \
    -out countRT.txt




