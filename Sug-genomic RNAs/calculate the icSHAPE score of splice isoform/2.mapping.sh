
D1=/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/5.map_human_genome/D1.Unmapped.out.mate1
D2=/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/5.map_human_genome/D2.Unmapped.out.mate1
N1=/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/5.map_human_genome/N1.Unmapped.out.mate1
N2=/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/5.map_human_genome/N2.Unmapped.out.mate1
T1=/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/5.map_human_genome/T1.Unmapped.out.mate1
T2=/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/5.map_human_genome/T2.Unmapped.out.mate1
OUT=/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/14.map_virus_splice-all

function mapping_virus(){
    input_fq=$1
    output_prefix=$2
    INDEX=/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/14.map_virus_splice-all/index/

    # Set --outSAMmultNmax -1 to output all alignments
    # Set scoreDelOpen=-99 to disable deletion
    # Set scoreInsOpen=-99 to disable insertion
    # Set scoreGap=-99 to disable junctions
    bsub -q Z-ZQF -e error -o log -n 20 \
        "STAR --runMode alignReads \
            --runThreadN 20 \
            --genomeDir ${INDEX} \
            --readFilesType Fastx \
            --readFilesIn ${input_fq} \
            --outFileNamePrefix ${output_prefix} \
            --outReadsUnmapped None \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMmode Full \
            --scoreDelOpen -99 \
            --scoreInsOpen -99 \
            --scoreGap -99 \
            --outSAMattributes All \
            --outSAMmultNmax -1 \
            --outFilterMultimapNmax -1 \
            --outFilterMismatchNmax 3 \
            --outFilterMultimapScoreRange 1 \
            --outFilterMismatchNoverLmax 0.2 \
            --outFilterIntronMotifs None \
            --outSJfilterReads All \
            --alignEndsType EndToEnd \
            --alignSoftClipAtReferenceEnds Yes \
            --chimOutType SeparateSAMold \
            --limitBAMsortRAM 9052078072"
}

mapping_virus ${D1} D1.
mapping_virus ${D2} D2.
mapping_virus ${N1} N1.
mapping_virus ${N2} N2.
mapping_virus ${T1} T1.
mapping_virus ${T2} T2.



