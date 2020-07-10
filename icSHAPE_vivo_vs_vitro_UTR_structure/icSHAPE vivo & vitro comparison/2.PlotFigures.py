
icSHAPE_rep1 = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/8.calcSHAPE/virus-w50-rep1.shape'
icSHAPE_rep2 = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/8.calcSHAPE/virus-w50-rep2.shape'
icSHAPE_vitro_rep1 = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/8.calcSHAPE/virus-vitro-w50-rep1.shape'
icSHAPE_vitro_rep2 = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/8.calcSHAPE/virus-vitro-w50-rep2.shape'

icSHAPE_rep1 = General.load_shape(icSHAPE_rep1)['NC_045512.2']
icSHAPE_rep2 = General.load_shape(icSHAPE_rep2)['NC_045512.2']
icSHAPE_vitro_rep1 = General.load_shape(icSHAPE_vitro_rep1)['NC_045512.2']
icSHAPE_vitro_rep2 = General.load_shape(icSHAPE_vitro_rep2)['NC_045512.2']

names = ['icSHAPE_invivo_rep1','icSHAPE_invivo_rep2','icSHAPE_invitro_rep1','icSHAPE_invitro_rep2']
df = General.init_pd_rect(4,4,names,names,1)

data = [ icSHAPE_rep1,icSHAPE_rep2,icSHAPE_vitro_rep1,icSHAPE_vitro_rep2 ]

for i in range(4):
    for j in range(i+1, 4):
        x = []
        y = []
        for s1,s2 in zip(data[i], data[j]):
            if 'NULL' not in (s1, s2):
                x.append(float(s1))
                y.append(float(s2))
        r, p = scipy.stats.pearsonr(x,y)
        df.iloc[i,j] = df.iloc[j,i] = r

sns.heatmap(df, cmap='YlGnBu', annot=True, fmt='.3f')
plt.tight_layout()
plt.savefig(join(HOME, "figs/icSHAPE-correlation.pdf"))
plt.close()

vivo1 = [float(x) for x in icSHAPE_rep1 if x!='NULL']
vivo2 = [float(x) for x in icSHAPE_rep2 if x!='NULL']
vitro1 = [float(x) for x in icSHAPE_vitro_rep1 if x!='NULL']
vitro2 = [float(x) for x in icSHAPE_vitro_rep2 if x!='NULL']

plt.boxplot([vivo1,vivo2,vitro1,vitro2])
plt.xticks([1,2,3,4],['icSHAPE invivo rep1', 'icSHAPE vivo rep2', 'icSHAPE invitro rep1', 'icSHAPE invitro rep2'])
plt.tight_layout()
plt.savefig(join(HOME, "figs/icSHAPE-vivo-vitro.pdf"))
plt.close()


print( np.mean(vivo1) )
print( np.mean(vivo2) )
print( np.mean(vitro1) )
print( np.mean(vitro2) )

ttest_ind = scipy.stats.ttest_ind
print( ttest_ind(vivo1+vivo2, vitro1+vitro2) )


