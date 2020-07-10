
#################
#### 这里主要是蛋白CDS区域的icSHAPE和翻译速率之间的关系
#################

importCommon()
import SARS2

sequence = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')['NC_045512.2']
vivo_shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-w50.shape')['NC_045512.2']
vitro_shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-vitro-w50.shape')['NC_045512.2']
vivo_shape = np.array([np.nan if d=='NULL' else float(d) for d in vivo_shape])
vitro_shape = np.array([np.nan if d=='NULL' else float(d) for d in vitro_shape])


############   6h
TE_6h = """
S   168.7152762 337.6610088 270.1319544
N   1517.045767 833.2091447 1033.033722
7a  14.86073794 45.28320747 151.0277983
""".strip().split('\n')
TE_6h = [ dataline.split() for dataline in TE_6h ]
columns = ['protein', 'rep1', 'rep2', 'rep3']
TE_6h = pd.DataFrame([ [dataline[0],float(dataline[1]),float(dataline[2]),float(dataline[3])] for dataline in TE_6h ], columns=columns)
TE_6h.index = TE_6h['protein']
TE_6h = TE_6h.drop('protein', axis=1)

############   10h
TE_10h = """
S       1117.576225 2265.556753 1082.508934
N       8230.723627 6489.96545  4374.156002
7a      111.5548113 388.0340874 760.7074596
""".strip().split('\n')
TE_10h = [ dataline.split() for dataline in TE_10h ]
columns = ['protein', 'rep1', 'rep2', 'rep3']
TE_10h = pd.DataFrame([ [dataline[0],float(dataline[1]),float(dataline[2]),float(dataline[3])] for dataline in TE_10h ], columns=columns)
TE_10h.index = TE_10h['protein']
TE_10h = TE_10h.drop('protein', axis=1)

############ 24h
TE_24h = """
S   1943.364451 3866.265217 1977.42077
N   15231.1338  11461.391   7149.18684
7a  183.7535429 568.1455106 1167.448313
""".strip().split('\n')
TE_24h = [ dataline.split() for dataline in TE_24h ]
columns = ['protein', 'rep1', 'rep2', 'rep3']
TE_24h = pd.DataFrame([ [dataline[0],float(dataline[1]),float(dataline[2]),float(dataline[3])] for dataline in TE_24h ], columns=columns)
TE_24h.index = TE_24h['protein']
TE_24h = TE_24h.drop('protein', axis=1)


########################
######################## 画图
########################

dataframe = General.init_pd_rect(5, 3, ['TE_6h','TE_10h','TE_24h','icSHAPE_vivo','icSHAPE_vitro'], ['S','N','7a'])
dataframe.loc['TE_6h'] = TE_6h.mean(axis=1)
dataframe.loc['TE_10h'] = TE_10h.mean(axis=1)
dataframe.loc['TE_24h'] = TE_24h.mean(axis=1)
dataframe.loc['icSHAPE_vivo'] = [
    np.nanmean(vivo_shape[21563:25384]), # S
    np.nanmean(vivo_shape[28274:29533]), # N
    np.nanmean(vivo_shape[27394:27759])  # 7a
]
dataframe.loc['icSHAPE_vitro'] = [
    np.nanmean(vitro_shape[21563:25384]), # S
    np.nanmean(vitro_shape[28274:29533]), # N
    np.nanmean(vitro_shape[27394:27759])  # 7a
]

dataframe = np.transpose(dataframe)
dataframe['protein'] = dataframe.index

plt.figure(figsize=(10, 6))

plt.subplot(2,3,1)
sns.scatterplot('icSHAPE_vivo', 'TE_6h', data=dataframe, hue='protein')
r, p = scipy.stats.pearsonr(dataframe['icSHAPE_vivo'], dataframe['TE_6h'])
plt.title(f"r={r:.3} p={p:.3}")

plt.subplot(2,3,2)
sns.scatterplot('icSHAPE_vivo', 'TE_10h', data=dataframe, hue='protein')
r, p = scipy.stats.pearsonr(dataframe['icSHAPE_vivo'], dataframe['TE_10h'])
plt.title(f"r={r:.3} p={p:.3}")

plt.subplot(2,3,3)
sns.scatterplot('icSHAPE_vivo', 'TE_24h', data=dataframe, hue='protein')
r, p = scipy.stats.pearsonr(dataframe['icSHAPE_vivo'], dataframe['TE_24h'])
plt.title(f"r={r:.3} p={p:.3}")

plt.subplot(2,3,4)
sns.scatterplot('icSHAPE_vitro', 'TE_6h', data=dataframe, hue='protein')
r, p = scipy.stats.pearsonr(dataframe['icSHAPE_vitro'], dataframe['TE_6h'])
plt.title(f"r={r:.3} p={p:.3}")

plt.subplot(2,3,5)
sns.scatterplot('icSHAPE_vitro', 'TE_10h', data=dataframe, hue='protein')
r, p = scipy.stats.pearsonr(dataframe['icSHAPE_vitro'], dataframe['TE_10h'])
plt.title(f"r={r:.3} p={p:.3}")

plt.subplot(2,3,6)
sns.scatterplot('icSHAPE_vitro', 'TE_24h', data=dataframe, hue='protein')
r, p = scipy.stats.pearsonr(dataframe['icSHAPE_vitro'], dataframe['TE_24h'])
plt.title(f"r={r:.3} p={p:.3}")

plt.tight_layout()
plt.savefig(join(HOME, 'figs/Structure & TE.pdf'))
plt.close()

########################
######################## 合起来画
########################


dataframe = dataframe.sort_values(by="TE_6h")
plt.plot(dataframe['icSHAPE_vivo'], dataframe['TE_6h'], '-', c=Colors.RGB['gray'])
sns.scatterplot('icSHAPE_vivo', 'TE_6h', data=dataframe, hue='protein')
plt.plot(dataframe['icSHAPE_vivo'], dataframe['TE_10h'], '-', c=Colors.RGB['gray'])
sns.scatterplot('icSHAPE_vivo', 'TE_10h', data=dataframe, hue='protein')
plt.plot(dataframe['icSHAPE_vivo'], dataframe['TE_24h'], '-', c=Colors.RGB['gray'])
sns.scatterplot('icSHAPE_vivo', 'TE_24h', data=dataframe, hue='protein')
data = dataframe.loc[:,('TE_6h','icSHAPE_vivo')].values.tolist() + dataframe.loc[:,('TE_10h','icSHAPE_vivo')].values.tolist() + dataframe.loc[:,('TE_24h','icSHAPE_vivo')].values.tolist()
r,p = scipy.stats.pearsonr([d[0] for d in data], [d[1] for d in data])
plt.title(f"r={r:.3} p={p}")
plt.savefig(join(HOME, 'figs/Structure & TE(invivo).pdf'))
plt.close()

dataframe = dataframe.sort_values(by="TE_6h")
plt.plot(dataframe['icSHAPE_vitro'], dataframe['TE_6h'], '-', c=Colors.RGB['gray'])
sns.scatterplot('icSHAPE_vitro', 'TE_6h', data=dataframe, hue='protein')
plt.plot(dataframe['icSHAPE_vitro'], dataframe['TE_10h'], '-', c=Colors.RGB['gray'])
sns.scatterplot('icSHAPE_vitro', 'TE_10h', data=dataframe, hue='protein')
plt.plot(dataframe['icSHAPE_vitro'], dataframe['TE_24h'], '-', c=Colors.RGB['gray'])
sns.scatterplot('icSHAPE_vitro', 'TE_24h', data=dataframe, hue='protein')
data = dataframe.loc[:,('TE_6h','icSHAPE_vitro')].values.tolist() + dataframe.loc[:,('TE_10h','icSHAPE_vitro')].values.tolist() + dataframe.loc[:,('TE_24h','icSHAPE_vivo')].values.tolist()
r,p = scipy.stats.pearsonr([d[0] for d in data], [d[1] for d in data])
plt.title(f"r={r:.3} p={p}")
plt.savefig(join(HOME, 'figs/Structure & TE(invitro).pdf'))
plt.close()



