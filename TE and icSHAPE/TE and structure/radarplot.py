
# Libraries
import matplotlib.pyplot as plt
import pandas as pd
from math import pi
import numpy as np

sequence = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')['NC_045512.2']
vivo_shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-w50.shape')['NC_045512.2']
vivo_shape = np.array([np.nan if d=='NULL' else float(d) for d in vivo_shape])

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

dataframe = General.init_pd_rect(5, 3, ['TE_6h','TE_10h','TE_24h','icSHAPE','splice_count'], ['S','N','7a'])
dataframe.loc['TE_6h'] = TE_6h.mean(axis=1)
dataframe.loc['TE_10h'] = TE_10h.mean(axis=1)
dataframe.loc['TE_24h'] = TE_24h.mean(axis=1)
dataframe.loc['icSHAPE'] = [
    np.nanmean(vivo_shape[21551:21562]), # S
    np.nanmean(vivo_shape[28254:28273]), # N
    np.nanmean(vivo_shape[27384:27393])  # 7a
]
dataframe.loc['splice_count'] = [ np.log10(1117576), np.log10(7012631), np.log10(688112) ]

#####  归一化

dataframe.loc['TE_6h'] /= 1000
dataframe.loc['TE_10h'] /= 1000
dataframe.loc['TE_24h'] /= 1000
dataframe.loc['icSHAPE'] *= 10

dataframe = np.transpose(dataframe)

# number of variable
categories=list(dataframe)
N = len(categories)

# Initialise the spider plot
ax = plt.subplot(111, polar=True)
for i in range(3):
    # We are going to plot the first line of the data frame.
    # But we need to repeat the first value to close the circular graph:
    values=dataframe.iloc[i].values.flatten().tolist()
    values += values[:1]
    
    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]
    
    # Draw one axe per variable + add labels labels yet
    plt.xticks(angles[:-1], categories, color='grey', size=8)
    
    # Draw ylabels
    ax.set_rlabel_position(0)
    plt.yticks([5,10,15,20], ["5","10","15","20"], color="grey", size=7)
    plt.ylim(0,16)
    
    # Plot data
    ax.plot(angles, values, linewidth=1, linestyle='solid')
    
    # Fill area
    ax.fill(angles, values, 'b', alpha=0.1)

plt.savefig(join(HOME, "figs/radar.pdf"))
plt.close()


