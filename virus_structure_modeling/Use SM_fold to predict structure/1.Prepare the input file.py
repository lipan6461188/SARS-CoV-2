
sequence = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')['NC_045512.2']
shape = General.load_shape('/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/8.calcSHAPE/virus-w50.shape')['NC_045512.2']
pdot = "."*len(sequence)

General.write_ct({'SARS2':sequence}, {'SARS2':pdot}, 'SARS2.ct')

print("\n".join([f'{i+1}\t-999' if d=='NULL' else str(i+1)+"\t"+str(d) for i,d in enumerate(shape)]), file=open('SARS2.shape','w'))

