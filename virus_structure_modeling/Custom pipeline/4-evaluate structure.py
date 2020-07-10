import General, Structure, Visual, SARS2
importCommon()

wuhan_id = 'NC_045512.2'
sequence = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')[wuhan_id]
shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-w50.shape')[wuhan_id]
dot = General.load_dot('/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/SARS2.dot')[wuhan_id][1]

####################
### 1. Overall AUC
####################

auc = General.calc_AUC_v2(dot, shape)
print(f"AUC = {auc:.3f}")

####################
### 2. icSHAPE distribution for single-stranded and double-stranded bases
####################

single_base_icSHAPE = [ float(s) for s,d in zip(shape, dot) if d=='.' and s!='NULL' ]
double_base_icSHAPE = [ float(s) for s,d in zip(shape, dot) if d!='.' and s!='NULL' ]

n1 = len(single_base_icSHAPE)
n2 = len(double_base_icSHAPE)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5), sharey=True)
ax.set_title(f'icSHAPE distribution\n n={n1},{n2}')
ax.set_ylabel('icSHAPE score')
data = [ single_base_icSHAPE, double_base_icSHAPE ]
colors = [ Colors.RGB['blue'], Colors.RGB['green'] ]
labels = [ 'Unpaired', 'Paired' ]
Figures.violinPlot(ax, data, labels, colors=colors)
fig.tight_layout()
plt.savefig(join(HOME, "figs/violin.pdf"))
plt.close()
#fig.show()

print(np.mean(single_base_icSHAPE))
print(np.mean(double_base_icSHAPE))
print(scipy.stats.mannwhitneyu(single_base_icSHAPE, double_base_icSHAPE, alternative='greater'))

# stems = Structure.find_stem(dot, max_stem_gap=1, min_stem_len=3)
# double_base_icSHAPE = []
# for ls,le,rs,re in stems:
#     double_base_icSHAPE += [ float(shape[i]) for i in range(ls, le-1) if shape[i]!='NULL' ]
#     double_base_icSHAPE += [ float(shape[i]) for i in range(rs, re-1) if shape[i]!='NULL' ]

# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5), sharey=True)
# ax.set_title('icSHAPE distribution')
# ax.set_ylabel('icSHAPE score')
# data = [ single_base_icSHAPE, double_base_icSHAPE ]
# colors = [ Colors.RGB['blue'], Colors.RGB['green'] ]
# labels = [ 'Unpaired', 'Paired' ]
# Figures.violinPlot(ax, data, labels, colors=colors)
# fig.tight_layout()
# fig.show()


####################
### 3. How many stem-loop structures
####################

stemloops = SARS2.find_stemloop(dot)
plt.boxplot([d[1]-d[0]+1 for d in stemloops])
n = len(stemloops)
plt.title(f"Stemloop length distribution (n={n})")
plt.show()


 