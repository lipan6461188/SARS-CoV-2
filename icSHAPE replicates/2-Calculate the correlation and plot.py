
importCommon()

def read_RT(infn):
    RT = {  
        'drt1': [0]*29903,
        'drt2': [0]*29903,
        'nrt1': [0]*29903,
        'nrt2': [0]*29903,
        'trt1': [0]*29903,
        'trt2': [0]*29903
    }
    IN = open(infn)
    for line in IN:
        if line[0]!='@':
            data = line.strip().split()
            if data[0]=='NC_045512.2' and data[1]=='+':
                pos = int(data[2])
                RT['drt1'][pos-1] = int(data[3])
                RT['drt2'][pos-1] = int(data[5])
                RT['nrt1'][pos-1] = int(data[7])
                RT['nrt2'][pos-1] = int(data[9])
                RT['trt1'][pos-1] = int(data[11])
                RT['trt2'][pos-1] = int(data[13])
    return RT

file = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/countRT.txt'
RT = read_RT(file)

labels = [ 'drt1','drt2','nrt1','nrt2','trt1','trt2' ]
cor = General.init_pd_rect(6,6,labels,labels,1)

for i in range(6):
    for j in range(1, 6):
        cor.iloc[i, j] = cor.iloc[j, i] = scipy.stats.pearsonr(RT[labels[i]], RT[labels[j]])[0]


sns.heatmap(cor, cmap='YlGnBu', annot=True, fmt=".3f")
plt.savefig(join(HOME, "figs/heatmap.pdf"))
plt.close()

plt.scatter(RT['drt1'], RT['drt2'])
plt.xlim(0, 20000)
plt.ylim(0, 20000)
plt.savefig(join(HOME, "figs/d_scatter.pdf"))
plt.close()

plt.scatter(RT['nrt1'], RT['nrt2'])
plt.xlim(0, 20000)
plt.ylim(0, 20000)
plt.savefig(join(HOME, "figs/n_scatter.pdf"))
plt.close()

plt.scatter(RT['trt1'], RT['trt2'])
plt.xlim(0, 20000)
plt.ylim(0, 20000)
plt.savefig(join(HOME, "figs/t_scatter.pdf"))
plt.close()






