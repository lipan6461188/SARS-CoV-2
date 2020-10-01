
workdir=/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing
V_SIZE=/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/index/SARS2_complement.size

icSHAPE-pipe calcSHAPE \
    -D ${workdir}/7.sam2tab/D1_virus.tab \
    -N ${workdir}/7.sam2tab/N1_virus.tab \
    -size ${V_SIZE} \
    -out ${workdir}/8.calcSHAPE/virus-w50-rep1.gTab \
    -noparam \
    -wsize 50 \
    -wstep 5 \
    -omc 0 \
    -mc 50

icSHAPE-pipe calcSHAPE \
    -D ${workdir}/7.sam2tab/D2_virus.tab \
    -N ${workdir}/7.sam2tab/N2_virus.tab \
    -size ${V_SIZE} \
    -out ${workdir}/8.calcSHAPE/virus-w50-rep2.gTab \
    -noparam \
    -wsize 50 \
    -wstep 5 \
    -omc 0 \
    -mc 50

icSHAPE-pipe genSHAPEToTransSHAPE \
    -i ${workdir}/8.calcSHAPE/virus-w50-rep1.gTab \
    -s ${V_SIZE} \
    -p 1 \
    -c 100 \
    -T 1 \
    -n 10 \
    -m 0 \
    -o ${workdir}/8.calcSHAPE/virus-w50-rep1.shape

icSHAPE-pipe genSHAPEToTransSHAPE \
    -i ${workdir}/8.calcSHAPE/virus-w50-rep2.gTab \
    -s ${V_SIZE} \
    -p 1 \
    -c 100 \
    -T 1 \
    -n 10 \
    -m 0 \
    -o ${workdir}/8.calcSHAPE/virus-w50-rep2.shape






icSHAPE-pipe calcSHAPE \
    -D ${workdir}/7.sam2tab/D1_virus.tab \
    -N ${workdir}/7.sam2tab/T1_virus.tab \
    -size ${V_SIZE} \
    -out ${workdir}/8.calcSHAPE/virus-vitro-w50-rep1.gTab \
    -noparam \
    -wsize 50 \
    -wstep 5 \
    -omc 0 \
    -mc 50

icSHAPE-pipe calcSHAPE \
    -D ${workdir}/7.sam2tab/D2_virus.tab \
    -N ${workdir}/7.sam2tab/T2_virus.tab \
    -size ${V_SIZE} \
    -out ${workdir}/8.calcSHAPE/virus-vitro-w50-rep2.gTab \
    -noparam \
    -wsize 50 \
    -wstep 5 \
    -omc 0 \
    -mc 50

icSHAPE-pipe genSHAPEToTransSHAPE \
    -i ${workdir}/8.calcSHAPE/virus-vitro-w50-rep1.gTab \
    -s ${V_SIZE} \
    -p 1 \
    -c 100 \
    -T 1 \
    -n 10 \
    -m 0 \
    -o ${workdir}/8.calcSHAPE/virus-vitro-w50-rep1.shape

icSHAPE-pipe genSHAPEToTransSHAPE \
    -i ${workdir}/8.calcSHAPE/virus-vitro-w50-rep2.gTab \
    -s ${V_SIZE} \
    -p 1 \
    -c 100 \
    -T 1 \
    -n 10 \
    -m 0 \
    -o ${workdir}/8.calcSHAPE/virus-vitro-w50-rep2.shape
