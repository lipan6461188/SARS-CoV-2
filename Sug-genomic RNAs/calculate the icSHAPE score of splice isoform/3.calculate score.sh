
cd /Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/14.map_virus_splice-all

bsub -q Z-ZQF -e error "icSHAPE-pipe sam2tab -in D1.Aligned.sortedByCoord.out.bam -out D1.tab"
bsub -q Z-ZQF -e error "icSHAPE-pipe sam2tab -in D2.Aligned.sortedByCoord.out.bam -out D2.tab"
bsub -q Z-ZQF -e error "icSHAPE-pipe sam2tab -in N1.Aligned.sortedByCoord.out.bam -out N1.tab"
bsub -q Z-ZQF -e error "icSHAPE-pipe sam2tab -in N2.Aligned.sortedByCoord.out.bam -out N2.tab"
bsub -q Z-ZQF -e error "icSHAPE-pipe sam2tab -in T1.Aligned.sortedByCoord.out.bam -out T1.tab"
bsub -q Z-ZQF -e error "icSHAPE-pipe sam2tab -in T2.Aligned.sortedByCoord.out.bam -out T2.tab"

bsub -q Z-ZQF -e error3 \
    "icSHAPE-pipe calcSHAPE \
        -D D1.tab,D2.tab \
        -N N1.tab,N2.tab \
        -size index/chrNameLength.txt \
        -wsize 50 \
        -out shape.gTab"

icSHAPE-pipe genSHAPEToTransSHAPE \
    -i shape.gTab \
    -s index/chrNameLength.txt \
    -c 100 \
    -o final.shape


cd only_cover_splice

awk '{split($1,arr,"_"); if($4>arr[1]){ print $0 }}' ../D1.tab > D1.tab
awk '{split($1,arr,"_"); if($4>arr[1]){ print $0 }}' ../D2.tab > D2.tab
awk '{split($1,arr,"_"); if($4>arr[1]){ print $0 }}' ../N1.tab > N1.tab
awk '{split($1,arr,"_"); if($4>arr[1]){ print $0 }}' ../N2.tab > N2.tab
awk '{split($1,arr,"_"); if($4>arr[1]){ print $0 }}' ../T1.tab > T1.tab
awk '{split($1,arr,"_"); if($4>arr[1]){ print $0 }}' ../T2.tab > T2.tab

bsub -q Z-ZQF -e error3 \
    "icSHAPE-pipe calcSHAPE \
        -D D1.tab,D2.tab \
        -N N1.tab,N2.tab \
        -size ../index/chrNameLength.txt \
        -wsize 50 \
        -out shape.gTab"

bsub -q Z-ZQF -e error3 \
    "icSHAPE-pipe calcSHAPE \
        -D D1.tab,D2.tab \
        -N T1.tab,T2.tab \
        -size ../index/chrNameLength.txt \
        -wsize 50 \
        -out shape_vitro.gTab"

icSHAPE-pipe genSHAPEToTransSHAPE \
    -i shape.gTab \
    -s ../index/chrNameLength.txt \
    -c 100 \
    -o final.shape

icSHAPE-pipe genSHAPEToTransSHAPE \
    -i shape_vitro.gTab \
    -s ../index/chrNameLength.txt \
    -c 100 \
    -o final_vitro.shape


