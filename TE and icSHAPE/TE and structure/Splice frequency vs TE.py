
importCommon()
import SARS2

sequence = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')['NC_045512.2']
vivo_shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-w50.shape')['NC_045512.2']
vitro_shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-vitro-w50.shape')['NC_045512.2']
vivo_shape = np.array([np.nan if d=='NULL' else float(d) for d in vivo_shape])
vitro_shape = np.array([np.nan if d=='NULL' else float(d) for d in vitro_shape])

splice_sites = """
64  28254   7012631 28273
65  21551   1117576 21562
66  27384   688112  27393
65  25380   618367  25392
65  27883   229749  27893
64  26467   215361  26522
69  26236   149744  26244
69  27040   38086   27201
68  28262   37740   28273
63  21549   17894   21562
63  28253   13079   28273
76  28266   9868    28273
60  27378   9320    27393
71  27761   8568    27824
61  28250   6457    28273
68  15776   6103    15811
70  27760   5982    69
69  27673   3866    27755
63  27881   3702    27893
61  28251   2946    28273
65  27383   2814    27393
66  27484   1901    27755
68  29153   1216    29221
76  21562   992 21562
74  21057   990 21070
76  27394   950 27755
83  21562   863 21562
74  21552   843 21562
68  29111   818 29221
67  24776   728 24809
70  22944   703 23003
70  22276   699 22358
68  28282   686 28283
82  26276   669 26522
75  28964   622 28972
67  26865   606 26895
70  22501   580 22613
71  27675   573 27755
62  28254   554 28273
69  24890   553 24890
70  29002   514 29078
68  27676   492 27755
71  26505   489 26522
68  21549   427 21562
75  26485   409 26522
63  25378   408 25392
70  21305   396 21316
69  18940   360 18940
68  27043   360 27201
78  25393   355 25404
69  29165   345 29221
68  27762   338 27824
72  21546   290 21562
64  26497   268 26522
63  5784    264 5858
69  23027   263 23045
68  25424   258 25456
77  26271   258 26522
69  19570   256 19630
75  27798   255 27824
57  21539   253 21562
70  20315   250 20329
63  26466   230 26522
68  26235   230 26244
64  27382   225 27393
67  2291    222 2455
70  14371   218 14380
70  2793    217 69
77  26480   212 26522
73  23027   210 23045
66  27298   204 27358
68  26289   188 26522
59  28249   186 28273
69  19111   185 19147
67  22404   185 22406
69  28164   184 28206
83  26250   183 26522
61  26464   179 26522
68  5788    178 5858
68  18050   175 18152
74  28265   175 28273
68  28994   158 29078
68  27866   153 27867
69  28206   151 28206
67  29378   149 29402
72  28570   146 28573
75  28267   144 28273
69  22558   144 22613
67  22502   143 22613
71  24892   137 24941
78  27754   136 27755
69  27758   135 27824
73  27777   134 27824
67  26925   134 27201
68  23364   133 23402
68  25490   133 25523
68  27473   133 27755
66  23849   131 24050
74  28891   128 28900
69  29590   126 29617
""".strip().split('\n')
splice_sites = [ d.split() for d in splice_sites ][:80]

def collect_SHAPE_Count(splice_sites, shape_array):
    func = np.nanmean
    df = []
    for p5,p3,count,start in splice_sites:
        p5,p3,count,start = int(p5),int(p3),int(count),int(start)
        e = max(p3, start)+30
        #if start>p3:
        #    e = start
        df.append( [func(shape_array[p3:e]), np.log10(count)] )
    return pd.DataFrame(df, columns=['score', 'count'])


vivo_df = collect_SHAPE_Count(splice_sites, vivo_shape)
#vitro_df = collect_SHAPE_Count(splice_sites, vitro_shape)
#VTD_df = collect_SHAPE_Count(splice_sites, vivo_shape-vitro_shape)

j = sns.jointplot(x='score', y='count', data=vivo_df, kind='kde', space=0, color=Colors.RGB['green']).annotate(scipy.stats.pearsonr).set_axis_labels("icSHAPE vivo", "log10(Count)")
plt.savefig(join(HOME, 'figs/SpliceCount_and_icSHAPE.pdf'))
plt.close()


###################################
#### 5'UTR区域的icSHAE值和Count的关系
###################################

###### 前10个、

splice_shape_fn = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/13.map_virus_splice/only_cover_splice/final.shape'
splice_shape = General.load_shape(splice_shape_fn)
df = []
for p5,p3,count,start in splice_sites[:10]:
    key = p5+"_"+p3
    shape_values = [float(it) for it in splice_shape[key][:int(p5)] if it!='NULL']
    if len(shape_values)<40: continue
    x = np.mean(shape_values)
    y =  np.log10(int(count))
    df.append([x, y, key])

df = pd.DataFrame(df, columns=['icSHAPE', 'log10(Count)', 'key'])
for i in range(df.shape[0]):
    x,y,label = df.iloc[i]
    plt.plot(x+0.02, y, '.', label=label)
    plt.text(x, y, label)

# 回归线
model = sklearn.linear_model.LinearRegression().fit(np.reshape(df['icSHAPE'].array, [-1,1]), np.reshape(df['log10(Count)'].array, [-1,1]))
k = model.coef_[0][0]
b = model.intercept_[0]
Xmin = 0.1# df.min(axis=0)[0]
Xmax = 0.25#df.max(axis=0)[0]
x = [ Xmin, Xmax ]
y = [ k*Xmin+b, k*Xmax+b ]
plt.plot(x, y, '-')

plt.xlim(0.1, 0.25)
r, p =scipy.stats.pearsonr(df['icSHAPE'], df['log10(Count)'])
plt.title(f"r={r:.3f}, p={p:.5f}")
plt.savefig(join(HOME, 'figs/SpliceCount_and_5picSHAPE.pdf'))
plt.close()

###### 前100个

splice_shape_fn = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/14.map_virus_splice-all/only_cover_splice/final.shape'
splice_shape = General.load_shape(splice_shape_fn)
df = []
for p5,p3,count,start in splice_sites:
    key = p5+"_"+p3
    shape_values = [float(it) for it in splice_shape[key][:int(p5)] if it!='NULL']
    if len(shape_values)<40: continue
    x = np.mean(shape_values)
    y =  np.log10(int(count))
    df.append([x, y, key])

df = pd.DataFrame(df, columns=['icSHAPE', 'log10(Count)', 'key'])

j = sns.jointplot(x='icSHAPE', y='log10(Count)', data=df, kind='kde', space=0, color=Colors.RGB['green']).annotate(scipy.stats.pearsonr).set_axis_labels("icSHAPE 5' vivo", "log10(Count)")
plt.savefig(join(HOME, 'figs/SpliceCount_and_5picSHAPE_all.pdf'))
plt.close()


j = sns.jointplot(x='icSHAPE', y='log10(Count)', data=df, kind='reg', space=0, color=Colors.RGB['green']).annotate(scipy.stats.pearsonr).set_axis_labels("icSHAPE 5' vivo", "log10(Count)")
plt.savefig(join(HOME, 'figs/SpliceCount_and_5picSHAPE_all.pdf'))
plt.close()

