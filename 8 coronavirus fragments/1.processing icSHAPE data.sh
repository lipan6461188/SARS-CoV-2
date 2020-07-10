
#########################################
### 由于发现孙磊把SARS2-C和SARS2-T搞反了所以需要重新mapping
#########################################


cd /Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/8UTR-icSHAPE/SARS2-C/OUT
samples=(D1 D2 N1 N2)
for sample in ${samples[@]};
do
    bsub -q Z-ZQF -n 20 -e error -o log \
        "icSHAPE-pipe cleanFq \
            -i 2.trim/${sample}.fastq \
            -o 3.rem_rRNA/${sample}.fastq \
            -x rRNA_index/SARS2-C \
            -p 20 \
            --mode Local \
            --sam 3.rem_rRNA/${sample}.sam"
done


cd /Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/8UTR-icSHAPE/SARS2-T-NL63_HKU1/OUT
samples=(D1 D2 N1 N2)
for sample in ${samples[@]};
do
    bsub -q Z-ZQF -n 20 -e error -o log \
        "icSHAPE-pipe cleanFq \
            -i 2.trim/${sample}.fastq \
            -o 3.rem_rRNA/${sample}.fastq \
            -x rRNA_index/SARS2-C_hNL63_hHKU1 \
            -p 20 \
            --mode Local \
            --sam 3.rem_rRNA/${sample}.sam"
done


#### 过滤掉SARS2-T中发生C突变的碱基

def filter(insam, outsam):
    IN = open(insam)
    OUT = open(outsam, 'w')
    for line in IN:
        if line[0]!='@':
            data = line.strip().split('\t')
            if data[2]=='SARS2-T' and int(data[3])<241 and int(data[3])+len(data[9])>241:
                attr22 = ""
                for attr in data[11:]:
                    if attr.startswith('MD:Z:'):
                        attr22 = attr
                assert attr22.startswith('MD:Z:'), attr22
                if 'C' in attr22:
                    continue
        print(line, file=OUT, end="")

filter('D1.sam', 'D1.sam2')
filter('D2.sam', 'D2.sam2')
filter('N1.sam', 'N1.sam2')
filter('N2.sam', 'N2.sam2')

mv D1.sam2 D1.sam
mv D2.sam2 D2.sam
mv N1.sam2 N1.sam
mv N2.sam2 N2.sam

#########################################
### 处理数据
#########################################


cd /Share/home/zhangqf8/sunlei/data/SARS2/20200623-inVitro-SRAS-huh7-8UTR/8UTR-icSHAPE

samples=(MERS-HKU9 SARS2-T-NL63_HKU1 SARS2-C SARS-HKU5)
for sample in ${samples[@]};
do
    awk 'substr($0,1,1)=="@"||$6!~/^[0-9]S/' ${sample}/OUT/3.rem_rRNA/D1.sam > ${sample}/OUT/3.rem_rRNA/D1-clean.sam
    awk 'substr($0,1,1)=="@"||$6!~/^[0-9]S/' ${sample}/OUT/3.rem_rRNA/D2.sam > ${sample}/OUT/3.rem_rRNA/D2-clean.sam
    awk 'substr($0,1,1)=="@"||$6!~/^[0-9]S/' ${sample}/OUT/3.rem_rRNA/N1.sam > ${sample}/OUT/3.rem_rRNA/N1-clean.sam
    awk 'substr($0,1,1)=="@"||$6!~/^[0-9]S/' ${sample}/OUT/3.rem_rRNA/N2.sam > ${sample}/OUT/3.rem_rRNA/N2-clean.sam
    
    icSHAPE-pipe sam2tab -in ${sample}/OUT/3.rem_rRNA/D1-clean.sam -out ${sample}/OUT/7.sam2tab/rRNA_D1-clean.tab
    icSHAPE-pipe sam2tab -in ${sample}/OUT/3.rem_rRNA/D2-clean.sam -out ${sample}/OUT/7.sam2tab/rRNA_D2-clean.tab
    icSHAPE-pipe sam2tab -in ${sample}/OUT/3.rem_rRNA/N1-clean.sam -out ${sample}/OUT/7.sam2tab/rRNA_N1-clean.tab
    icSHAPE-pipe sam2tab -in ${sample}/OUT/3.rem_rRNA/N2-clean.sam -out ${sample}/OUT/7.sam2tab/rRNA_N2-clean.tab

    icSHAPE-pipe calcSHAPE \
        -D ${sample}/OUT/7.sam2tab/rRNA_D1-clean.tab,${sample}/OUT/7.sam2tab/rRNA_D2-clean.tab \
        -N ${sample}/OUT/7.sam2tab/rRNA_N1-clean.tab,${sample}/OUT/7.sam2tab/rRNA_N2-clean.tab \
        -size ${sample}/OUT/rRNA_index/rRNA.len \
        -noparam \
        -wsize 50 \
        -out ${sample}/OUT/8.calcGenomeSHAPE/virus.gTab

    icSHAPE-pipe genSHAPEToTransSHAPE \
        -i ${sample}/OUT/8.calcGenomeSHAPE/virus.gTab \
        -s ${sample}/OUT/rRNA_index/rRNA.len \
        -o ${sample}/OUT/9.transSHAPE/virus.shape
done


