
importCommon()

def read_UTR_icSHAPE():
    root = "/Share/home/zhangqf8/sunlei/data/SARS2/20200623-inVitro-SRAS-huh7-8UTR/8UTR-icSHAPE"
    icSHAPE = {}
    icSHAPE.update( General.load_shape(join(root, "MERS-HKU9/OUT/9.transSHAPE/virus.shape")) )
    icSHAPE.update( General.load_shape(join(root, "SARS2-T-NL63_HKU1/OUT/9.transSHAPE/virus.shape")) )
    icSHAPE.update( General.load_shape(join(root, "SARS2-C/OUT/9.transSHAPE/virus.shape")) )
    icSHAPE.update( General.load_shape(join(root, "SARS-HKU5/OUT/9.transSHAPE/virus.shape")) )
    return icSHAPE

utr_icSHAPE = read_UTR_icSHAPE()
utr_sequences = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/sequences/utr_sequences.fa')

utr5_frag = {
    'SARS': 395,
    'SARS2-C': 395,
    'SARS2-T': 395,
    'bHKU5': 445,
    'MERS': 460,
    'hNL63': 330,
    'hHKU1': 330,
    'bHKU9': 330
}

utr3_frag = {
    'SARS': 330,
    'SARS2-C': 330,
    'SARS2-T': 330,
    'bHKU5': 290,
    'MERS': 290,
    'hNL63': 300,
    'hHKU1': 200,
    'bHKU9': 265
}

groups = [
    ['SARS', 'SARS2-C', 'SARS2-T'],
    ['bHKU5', 'MERS'],
    ['hNL63'],
    ['hHKU1'],
    ['bHKU9']
]

###### ticks is base

plt.figure(figsize=(12,8))
ii = 0
for group in groups:
    aligned_5utrs = Structure.multi_alignment([ utr_sequences[ID][:utr5_frag[ID]] for ID in group ])
    raw_shape = [ utr_icSHAPE[ID][:utr5_frag[ID]] for ID in group ]
    aligned_shape = [ Structure.shape_to_alignSHAPE(shape, aligned_5utrs[i]) for i,shape in enumerate(raw_shape) ]
    for i in range(len(group)):
        ii += 1
        plt.subplot(8,1,ii)
        values = [ 1 if v=='NULL' else float(v) for v in aligned_shape[i] ]
        colors = [ Colors.RGB['gray'] if v=='NULL' else Colors.RGB['blue'] for v in aligned_shape[i] ]
        plt.bar(range(len(aligned_5utrs[i])), values, color=colors)
        plt.xticks(range(len(aligned_5utrs[i])), list(aligned_5utrs[i]), fontsize=3)
        plt.title(group[i])

plt.tight_layout()
plt.savefig(join(HOME, 'figs/utr_profile.pdf'))
plt.close()

###### ticks is number

bar_colors = {
    'SARS2-C': '#f29600',
    'SARS2-T': '#ed7648',
    'SARS': '#dcb656',
    'bHKU5': '#2bb17f',
    'MERS': '#6cb94f',
    'hNL63': '#c8bb9b',
    'hHKU1': '#5ea4be',
    'bHKU9': '#5e7bbb'
}

######  5' UTR

plt.figure(figsize=(12,8))
ii = 0
for group in groups:
    aligned_5utrs = Structure.multi_alignment([ utr_sequences[ID][:utr5_frag[ID]] for ID in group ])
    raw_shape = [ utr_icSHAPE[ID][:utr5_frag[ID]] for ID in group ]
    aligned_shape = [ Structure.shape_to_alignSHAPE(shape, aligned_5utrs[i]) for i,shape in enumerate(raw_shape) ]
    for i in range(len(group)):
        ii += 1
        plt.subplot(8,1,ii)
        values = [ 1 if v=='NULL' else float(v) for v in aligned_shape[i] ]
        colors = [ Colors.RGB['gray'] if v=='NULL' else bar_colors[group[i]] for v in aligned_shape[i] ]
        plt.bar(range(len(aligned_5utrs[i])), values, color=colors)
        x = []; f = []
        j = 0
        for k,b in enumerate(aligned_5utrs[i]):
            if b!='-':
                j += 1
                if j%25==0:
                    x.append(k)
                    f.append(str(j))
        plt.xticks(x, f, fontsize=7)
        plt.title(group[i])

plt.tight_layout()
plt.savefig(join(HOME, 'figs/utr_profile_5p.pdf'))
plt.close()

######  3' UTR

plt.figure(figsize=(12,8))
ii = 0
for group in groups:
    aligned_3utrs = Structure.multi_alignment([ utr_sequences[ID][-utr3_frag[ID]:] for ID in group ])
    raw_shape = [ utr_icSHAPE[ID][-utr3_frag[ID]:] for ID in group ]
    aligned_shape = [ Structure.shape_to_alignSHAPE(shape, aligned_3utrs[i]) for i,shape in enumerate(raw_shape) ]
    for i in range(len(group)):
        ii += 1
        plt.subplot(8,1,ii)
        values = [ 1 if v=='NULL' else float(v) for v in aligned_shape[i] ]
        colors = [ Colors.RGB['gray'] if v=='NULL' else bar_colors[group[i]] for v in aligned_shape[i] ]
        plt.bar(range(len(aligned_3utrs[i])), values, color=colors)
        x = []; f = []
        j = 0
        for k,b in enumerate(aligned_3utrs[i]):
            if b!='-':
                j += 1
                if j%25==0:
                    x.append(k)
                    f.append(str(j))
        plt.xticks(x, f, fontsize=7)
        plt.title(group[i])

plt.tight_layout()
plt.savefig(join(HOME, 'figs/utr_profile_3p.pdf'))
plt.close()

# ##########################
# ####  画pair-wise sequence similarity
# ##########################

def seq_similarity(seq1, seq2):
    alignments = Structure.multi_alignment([seq1, seq2])
    sample_count = 0
    total = len(alignments[0])
    for b1,b2 in zip(alignments[0],alignments[1]):
        if b1==b2 and b1!='-':
            sample_count += 1
    return sample_count/total

def shape_similarity(seq1, seq2, shape1, shape2):
    alignments = Structure.multi_alignment([seq1, seq2])
    aligned_shape1 = Structure.shape_to_alignSHAPE(shape1, alignments[0])
    aligned_shape2 = Structure.shape_to_alignSHAPE(shape2, alignments[1])
    x = []
    y = []
    for s1,s2 in zip(aligned_shape1, aligned_shape2):
        if s1!='NULL' and s2!='NULL':
            x.append( float(s1) )
            y.append( float(s2) )
    r, p = scipy.stats.pearsonr(x, y)
    return r


#### Sequence similarity

# 5' UTR
UTR5_seq = { key:utr_sequences[key][:utr5_frag[key]] for key in utr5_frag }
keys = ['SARS', 'SARS2-C', 'SARS2-T', 'bHKU5', 'MERS', 'hNL63', 'hHKU1', 'bHKU9' ]
coor = General.init_pd_rect(8, 8, keys, keys, 1)

for i in range(8):
    for j in range(i+1, 8):
        coor.iloc[i, j] = coor.iloc[j, i] = seq_similarity(UTR5_seq[keys[i]], UTR5_seq[keys[j]])

sns.heatmap(coor, cmap="YlGnBu",  annot=True, fmt='.3f')
plt.tight_layout()
plt.savefig(join(HOME, 'figs/sequence_similarity_5p.pdf'))
plt.close()

# 3' UTR
UTR3_seq = { key:utr_sequences[key][-utr3_frag[key]:] for key in utr3_frag }
keys = ['SARS', 'SARS2-C', 'SARS2-T', 'bHKU5', 'MERS', 'hNL63', 'hHKU1', 'bHKU9' ]
coor = General.init_pd_rect(8, 8, keys, keys, 1)

for i in range(8):
    for j in range(i+1, 8):
        coor.iloc[i, j] = coor.iloc[j, i] = similarity(UTR3_seq[keys[i]], UTR3_seq[keys[j]])

sns.heatmap(coor, cmap="YlGnBu",  annot=True, fmt='.3f')
plt.tight_layout()
plt.savefig(join(HOME, 'figs/sequence_similarity_3p.pdf'))
plt.close()

#### Shape similarity

# 5' UTR
UTR5_seq = { key:utr_sequences[key][:utr5_frag[key]] for key in utr5_frag }
UTR5_shape = { key:utr_icSHAPE[key][:utr5_frag[key]] for key in utr5_frag }
keys = ['SARS', 'SARS2-C', 'SARS2-T', 'bHKU5', 'MERS', 'hNL63', 'hHKU1', 'bHKU9' ]
coor = General.init_pd_rect(8, 8, keys, keys, 1)

for i in range(8):
    for j in range(i+1, 8):
        coor.iloc[i, j] = coor.iloc[j, i] = shape_similarity(UTR5_seq[keys[i]], UTR5_seq[keys[j]], UTR5_shape[keys[i]], UTR5_shape[keys[j]])

sns.heatmap(coor, cmap="YlGnBu",  annot=True, fmt='.3f')
plt.tight_layout()
plt.savefig(join(HOME, 'figs/shape_similarity_5p.pdf'))
plt.close()

# 3' UTR
UTR3_seq = { key:utr_sequences[key][-utr3_frag[key]:] for key in utr3_frag }
UTR3_shape = { key:utr_icSHAPE[key][-utr3_frag[key]:] for key in utr3_frag }
keys = ['SARS', 'SARS2-C', 'SARS2-T', 'bHKU5', 'MERS', 'hNL63', 'hHKU1', 'bHKU9' ]
coor = General.init_pd_rect(8, 8, keys, keys, 1)

for i in range(8):
    for j in range(i+1, 8):
        coor.iloc[i, j] = coor.iloc[j, i] = shape_similarity(UTR3_seq[keys[i]], UTR3_seq[keys[j]], UTR3_shape[keys[i]], UTR3_shape[keys[j]])

sns.heatmap(coor, cmap="YlGnBu",  annot=True, fmt='.3f')
plt.tight_layout()
plt.savefig(join(HOME, 'figs/shape_similarity_3p.pdf'))
plt.close()





