
import General, Structure, Visual, SARS2
importCommon()

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
####  比较invivo 和 invitro之间差异的显著性
####################################

sites_dataframe = []
for i in range(len(icSHAPE_rep1)):
    a1,a2,b1,b2 = icSHAPE_rep1[i],icSHAPE_rep2[i],icSHAPE_vitro_rep1[i],icSHAPE_vitro_rep2[i]
    a1 = float(a1) if a1!='NULL' else np.nan
    a2 = float(a2) if a2!='NULL' else np.nan
    b1 = float(b1) if b1!='NULL' else np.nan
    b2 = float(b2) if b2!='NULL' else np.nan
    if not np.any(np.isnan(np.array([a1,a2,b1,b2]))):
        s1,p1 = scipy.stats.ttest_ind([a1,a2], [b1,b2])
        if np.isnan(p1):
            sites_dataframe.append([a1,a2,b1,b2,1.0])
        else:
            sites_dataframe.append([a1,a2,b1,b2,round(p1,8)])
    else:
        sites_dataframe.append([a1,a2,b1,b2,np.nan])

sites_dataframe = pd.DataFrame(sites_dataframe, columns=['invivo_rep1', 'invivo_rep2', 'invitro_rep1', 'invitro_rep2', 'Pvalue'])
sites_dataframe['invivo_mean'] = sites_dataframe.loc[:, ('invivo_rep1', 'invivo_rep2')].mean(axis=1)
sites_dataframe['invitro_mean'] = sites_dataframe.loc[:, ('invitro_rep1', 'invitro_rep2')].mean(axis=1)
sites_dataframe['VTD'] = sites_dataframe['invivo_mean'] - sites_dataframe['invitro_mean']
sites_dataframe = sites_dataframe.drop('invivo_mean', axis=1).drop('invitro_mean', axis=1)

pvalue_sites = ~np.isnan(sites_dataframe['Pvalue'].values)
pvals = sites_dataframe['Pvalue'].values[pvalue_sites]
pvals_adj = statsmodels.stats.multitest.multipletests(pvals, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
sites_dataframe['Pvalue_adj'] = np.nan
sites_dataframe.loc[pvalue_sites, 'Pvalue_adj'] = pvals_adj[1]

sites_dataframe.insert(0, "Pos", range(1, sites_dataframe.shape[0]+1))
sites_dataframe.insert(1, "Nuc", list(seq))

###############################
###  寻找显著的区域
###############################

def search_sig_regions(p_value_list, abs_diff_list, min_sig_num=3, min_region_win=5):
    window = []
    for i in range(len(p_value_list)-min_region_win):
        j = i+min_region_win
        p_list = p_value_list[i:j]
        d_list = abs_diff_list[i:j]
        if np.sum(np.isnan(p_list))==0:
            condition1 = (p_list<=0.05) & (d_list>=0.2)
            if np.sum(condition1)>=min_sig_num and condition1[0]==True and condition1[-1]==True:
                window.append([i, j])
    i = 0
    while i<len(window)-1:
        if window[i][1]>=window[i+1][0]:
            window[i] = [window[i][0], max(window[i][1], window[i+1][1]) ]
            del window[i+1]
        else:
            i += 1
    return window

window = search_sig_regions(sites_dataframe.Pvalue_adj.values, np.abs(sites_dataframe['VTD']).values)

for ws,we in window[:10]+window[-10:]:
    print( Colors.f(f"{ws}-{we}", bc='yellow') )
    print(sites_dataframe.iloc[ws:we, :])

OUT = open(join(HOME, 'figs/vivo-vitro-diff.csv'), 'w')
print("ID,Start,End,Annot,Shape vivo rep1,Shape vivo rep2,Shape vitro rep1,Shape vitro rep2", file=OUT)
for ws,we in window:
    shape_vivo_rep1 = " ".join(icSHAPE_rep1[ws:we])
    shape_vivo_rep2 = " ".join(icSHAPE_rep2[ws:we])
    shape_vitro_rep1 = " ".join(icSHAPE_vitro_rep1[ws:we])
    shape_vitro_rep2 = " ".join(icSHAPE_vitro_rep2[ws:we])
    annot = SARS2.annotate_region(ws+1, we)
    print(f"NC_045512.2,{ws+1},{we},{annot},{shape_vivo_rep1},{shape_vivo_rep2},{shape_vitro_rep1},{shape_vitro_rep2}", file=OUT)

OUT.close()


########################## 统计每一个ORF中差异区域的比例


diff_bases_VTD_gt0 = {}
diff_bases_VTD_lt0 = {}
for ws,we in window:
    annots = SARS2.annotate_region(ws+1, we).split()
    diff = [float(a)+float(b)-float(c)-float(d) for a,b,c,d in zip(icSHAPE_rep1[ws:we], icSHAPE_rep2[ws:we],icSHAPE_vitro_rep1[ws:we],icSHAPE_vitro_rep2[ws:we])]
    diff = np.mean(diff)
    for annot in annots:
        if diff>0:
            diff_bases_VTD_gt0[annot] = diff_bases_VTD_gt0.get(annot, 0) + (we-ws)
        else:
            diff_bases_VTD_lt0[annot] = diff_bases_VTD_lt0.get(annot, 0) + (we-ws)

ORF = SARS2.read_ORF()

df_diff = []
for name in ORF:
    length = ORF[name][1] - ORF[name][0] + 1
    diff1 = round(diff_bases_VTD_gt0.get(name,0)/length, 3)
    diff2 = round(diff_bases_VTD_lt0.get(name,0)/length, 3)
    df_diff.append([name, diff1, diff2])

df_diff = pd.DataFrame(df_diff, columns=['Region', 'invivo>invitro', 'invivo<invitro'])
df_diff.to_csv(join(HOME, "figs/diff_region_annotation.csv"), index=False)


###############################
###  在画差异窗口的同时把显著性窗口也标记上去
###############################

highlight_region = [0]*len(seq)
for start,end in window:
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





