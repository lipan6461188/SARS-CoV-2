
import General, Structure, Visual, SARS2
importCommon()
from scipy.stats import binom

wuhan_id = 'NC_045512.2'
sequence = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')[wuhan_id]
vivoshape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-w50.shape')[wuhan_id]
vitroshape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-vitro-w50.shape')[wuhan_id]

icSHAPE_rep1 = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/8.calcSHAPE/virus-w50-rep1.shape'
icSHAPE_rep2 = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/8.calcSHAPE/virus-w50-rep2.shape'
icSHAPE_vitro_rep1 = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/8.calcSHAPE/virus-vitro-w50-rep1.shape'
icSHAPE_vitro_rep2 = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/8.calcSHAPE/virus-vitro-w50-rep2.shape'

icSHAPE_rep1 = General.load_shape(icSHAPE_rep1)['NC_045512.2']
icSHAPE_rep2 = General.load_shape(icSHAPE_rep2)['NC_045512.2']
icSHAPE_vitro_rep1 = General.load_shape(icSHAPE_vitro_rep1)['NC_045512.2']
icSHAPE_vitro_rep2 = General.load_shape(icSHAPE_vitro_rep2)['NC_045512.2']

wuhan_id = 'NC_045512.2'
seq = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')[wuhan_id]

####################################
####  先转成numpy的array
####################################

vivoshape = np.array([ np.nan if d=='NULL' else float(d) for d in vivoshape ])
vitroshape = np.array([ np.nan if d=='NULL' else float(d) for d in vitroshape ])
VTD = np.abs(vivoshape - vitroshape)

vivoshape_rep1 = np.array([ np.nan if d=='NULL' else float(d) for d in icSHAPE_rep1 ])
vivoshape_rep2 = np.array([ np.nan if d=='NULL' else float(d) for d in icSHAPE_rep2 ])
vivoshape_vitro_rep1 = np.array([ np.nan if d=='NULL' else float(d) for d in icSHAPE_vitro_rep1 ])
vivoshape_vitro_rep2 = np.array([ np.nan if d=='NULL' else float(d) for d in icSHAPE_vitro_rep2 ])

vivo_variation = np.abs(vivoshape_rep1 - vivoshape_rep2)
vitro_variation = np.abs(vivoshape_vitro_rep1 - vivoshape_vitro_rep2)

ingroup_variation = (vivo_variation+vitro_variation) / 2


# 画图

cond = ~ (np.isnan(ingroup_variation) | np.isnan(VTD))
sns.scatterplot(ingroup_variation[cond], VTD[cond])
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel("Within invivo or invitro")
plt.ylabel("Between invivo and invitro")
plt.savefig(join(HOME, 'figs/invivo_invitro_scatter.pdf'))
plt.show()


####################################
####  求所需要的差异位点个数
####################################

ingroup_datalist = vivo_variation[~np.isnan(vivo_variation)].tolist() + vitro_variation[~np.isnan(vitro_variation)].tolist()
quantile = np.quantile(ingroup_datalist, 0.95)

window_size = 10
print(1-binom.cdf(2, window_size, 0.05)) # 小于等于3个差异碱基出现的概率


####################################
####  Call windows
####################################

raw_windows = []
for i in range(0, len(VTD)-window_size):
    small_window = VTD[i:i+window_size]
    if sum(small_window>quantile)>=2:
        raw_windows.append([i, i+window_size])

print(len(raw_windows))

# 按照窗口差异值排序，取前10%的

diff_windows = []
for i in range(len(raw_windows)):
    s,e = raw_windows[i]
    if sum(~np.isnan(VTD[s:e])) >= 5:
        diff = np.nanmean(VTD[s:e])
        diff_windows.append( [s,e,round(diff, 3)] )

sorted_diff_windows = sorted(diff_windows, key=lambda x: x[2])
top_num = int(len(sorted_diff_windows) * 0.1)
windows = sorted(sorted_diff_windows[-top_num:], key=lambda x: x[0])

print(len(windows))

# 合并窗口

i = 1
while i<len(windows):
    if windows[i-1][1]>=windows[i][0]:
        windows[i-1][1] = windows[i][1]
        del windows[i]
    else:
        i += 1

print(len(windows))
print(windows)

####################################
#### 输出保存这些区域
####################################

OUT = open("/tmp/diff.csv", 'w')
print("ID,Start,End,ORF,Binomial P-value,L1 distance", file=OUT)
for s, e, diff in windows:
    test = e - s
    success = sum(VTD[s:e]>quantile)
    pvalue = 1-binom.cdf(success, test, 0.05)
    L1 = np.nanmean(VTD[s:e])
    region = SARS2.annotate_region(s, e)
    print(f"NC_045512.2,{s+1},{e},{region},{pvalue:.5e},{L1:.3f}", file=OUT)
    if pvalue>0.05:
        print("warning")

OUT.close()

####################################
#### 注释这些区域
####################################

diff_r = {}
for s,e,_ in windows:
    annot = SARS2.annotate_region(s, e)
    for a in annot.split():
        try:
            diff_r[a]
        except:
            diff_r[a] = 0
        diff_r[a] = diff_r.get(a, 0) + (e-s)

df = []
for name in SARS2.read_ORF():
    length = SARS2.read_ORF()[name][1]-SARS2.read_ORF()[name][0]
    df.append( [ name, diff_r[name]/length ] )

df.sort(key=lambda x: x[1])

plt.bar(range(len(df)), [d[1] for d in df])
plt.xticks(range(len(df)), [d[0] for d in df], rotation=90)
plt.ylabel("Diff ratio")
plt.tight_layout()
plt.show()


####################################
####  和原来的结果做overlap
####################################

df = pd.read_csv("/Share2/home/zhangqf7/figs/vivo-vitro-diff.csv")
pan_pipe_diff_windows = df.loc[:, ('Start', 'End')].values

common = 0
for w in windows:
    for b in pan_pipe_diff_windows:
        if w[0]<b[1]-5 and b[0]<w[1]-5:
            common += 1
            break

print(common)

set1 = set(range(len(windows)))
set2 = set(range(len(windows)-common, len(windows)-common+len(pan_pipe_diff_windows)))
venn2( [set1, set2], set_labels=['WenZe', 'Pan'] )
plt.savefig(join(HOME, 'figs/method-overlap.pdf'))
plt.show()


###############################
###  画barplot
###############################

highlight_region = [0]*len(seq)
for start,end,_ in windows:
    for i in range(start, end):
        highlight_region[i] = 1


VTD = []
colors = []
for v1,v2 in zip(vivoshape, vitroshape):
    if v1!='NULL' and v2!='NULL':
        diff = float(v1) - float(v2)
        VTD.append(diff)
        if diff<0:
            colors.append(Colors.RGB['blue'])
        else:
            colors.append(Colors.RGB['green'])
    else:
        VTD.append(1)
        colors.append(Colors.RGB['gray'])

plt.subplot(2,1,1)
plt.bar(range(300), VTD[:300], color=colors[:300])
plt.subplot(2,1,2)
plt.bar(range(300), highlight_region[:300], color=Colors.RGB['gray'])
plt.savefig(join(HOME, "figs/vtd-5utr.pdf"))
plt.show()

plt.subplot(2,1,1)
plt.bar(range(29535,29903), VTD[29535:], color=colors[29535:])
plt.subplot(2,1,2)
plt.bar(range(29535,29903), highlight_region[29535:], color=Colors.RGB['gray'])
plt.savefig(join(HOME, "figs/vtd-3utr.pdf"))
plt.show()


plt.subplot(2,1,1)
plt.bar(range(1265,1507), VTD[1264:1506], color=colors[1264:1506])
plt.subplot(2,1,2)
plt.bar(range(1265,1507), highlight_region[1264:1506], color=Colors.RGB['gray'])
plt.savefig(join(HOME, "figs/vtd-cds.pdf"))
plt.show()




