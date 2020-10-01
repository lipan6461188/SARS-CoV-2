 
importCommon()
import SARS2

sequence = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')['NC_045512.2']
vivo_shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-w50.shape')['NC_045512.2']
vivo_shape = np.array([np.nan if d=='NULL' else float(d) for d in vivo_shape])

##### Top 100 RNAs
top100_fn = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/14.map_virus_splice-all/only_cover_splice/final.shape'
top100_shape = General.load_shape(top100_fn)
top100_shape = { k:np.array([np.nan if d=='NULL' else float(d) for d in top100_shape[k]]) for k in top100_shape }

top100_fa = General.load_fasta('/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/14.map_virus_splice-all/index/splice_sequence.fa')

cell_splice = pd.read_csv("/Share2/home/zhangqf7/lipan/SARS2/Revision/Cell-RNA-Seq/Cell-splice.csv")
media9 = pd.read_csv("/Share2/home/zhangqf7/lipan/SARS2/Revision/TE-Structure/TE_data.csv")

################ Read all isoforms data

class Isoform:
    def __init__(self, p5, p3, translation_start, count, annot):
        self.name = f"{p5}_{p3}"
        self.p5 = int(p5)
        self.p3 = int(p3)
        self.tstart = int(translation_start)
        self.count = int(count)
        self.annot = annot if type(annot)==str else 'NULL'
    
    def __repr__(self):
        return f"{self.p5}-{self.p3} start:{self.tstart} count:{self.count} annot:{self.annot}"
    
    def __str__(self):
        return f"{self.p5}-{self.p3} start:{self.tstart} count:{self.count} annot:{self.annot}"

def build_isoform(cell_splice, head=80):
    names = cell_splice.index
    columns = [ "5' site", "3' site", "startpos", "count", "startmatch_name"]
    values = cell_splice.loc[:, columns].values.tolist()
    Isoform_list = []
    i = 0
    for data in values[:head]:
        obj = Isoform(*data)
        Isoform_list.append(obj)
        i += 1
    return Isoform_list

top80_isoform_list = build_isoform(cell_splice, head=80)

def get_canonical_isoform(Isoform_list):
    annotate_sites = """
    266 ORF1ab
    21563 S
    25393 ORF3a
    26245 E
    26523 M
    27202 ORF6
    27394 ORF7a
    27756 ORF7b
    27894 ORF8
    28274 N
    29558 ORF10
    """.strip().split('\n')
    # 27756 ORF7b
    annotate_sites = [ d.split() for d in annotate_sites ]
    annotate_sites = { d[1]:int(d[0]) for d in annotate_sites }
    
    canonical_isoform = []
    proteins = set()
    for isoform_obj in Isoform_list:
        annot = isoform_obj.annot
        if annot not in ("NULL", 'pp1ab') and annot not in proteins:
            if isoform_obj.tstart+1 == annotate_sites[annot]:
                proteins.add(annot)
                canonical_isoform.append(isoform_obj)
    
    print(len(canonical_isoform))
    return canonical_isoform

canonical_isoform_list = get_canonical_isoform(top80_isoform_list)
# 去除ORF7b
canonical_isoform_list = [ d for d in canonical_isoform_list if d.annot!='ORF7b' ]
print(len(canonical_isoform_list))

def get_region_icSHAPE_mean_relative_Splice(input_Isoform_list, start, length, TE=None):
    x = []
    y = []
    names = []
    for obj in input_Isoform_list:
        if obj.name in top100_shape:
            if obj.annot == 'NULL':
                continue
            if TE and obj.name not in TE:
                continue
            p5, p3 = obj.name.split('_')
            p5, p3 = int(p5), int(p3)
            local_shape = top100_shape[obj.name][p5+start:p5+start+length]
            local_shape = local_shape[ ~np.isnan(local_shape) ]
            if len(local_shape)>=10:
                mean_shape = np.mean(local_shape)
                x.append(mean_shape)
                if TE is None:
                    log10count = np.log10(obj.count)
                else:
                    log10count = TE[obj.name]
                y.append(log10count)
                names.append(obj.annot)
    return x, y, names

def get_region_icSHAPE_mean_relative_Translation(input_Isoform_list, start, length, TE=None):
    x = []
    y = []
    names = []
    for obj in input_Isoform_list:
        if obj.name in top100_shape:
            if obj.annot == 'NULL':
                continue
            if TE and obj.name not in TE:
                continue
            p5, p3 = obj.name.split('_')
            p5, p3 = int(p5), int(p3)
            if obj.tstart<=obj.p5:
                tstart = obj.tstart
            else:
                assert obj.tstart>=obj.p3
                tstart = obj.tstart - obj.p3 + 1 + obj.p5
            local_shape = top100_shape[obj.name][tstart+start:tstart+start+length]
            local_shape = local_shape[ ~np.isnan(local_shape) ]
            if len(local_shape)>=10:
                mean_shape = np.mean(local_shape)
                log10count = np.log10(obj.count)
                x.append(mean_shape)
                if TE is None:
                    log10count = np.log10(obj.count)
                else:
                    log10count = TE[obj.name]
                y.append(log10count)
                names.append(obj.annot)
    return x, y, names

def calc_count_shape_corr(input_Isoform_list, func=None, TE=None):
    start_list = [-50, -40, -30, -20, -10, 0, 10, 20]
    length_list = [10, 20, 30, 40, 50, 60, 70]
    P_value_df = General.init_pd_rect(len(start_list), len(length_list), start_list, length_list, np.nan)
    Corr_df = General.init_pd_rect(len(start_list), len(length_list), start_list, length_list, np.nan)
    Datasize_df = General.init_pd_rect(len(start_list), len(length_list), start_list, length_list, 0)
    for start in start_list:
        for length in length_list:
            x, y, names = func(input_Isoform_list, start, length, TE=TE)
            if len(x)>5:
                corr, pvalue = scipy.stats.spearmanr(x, y)
                P_value_df.loc[start, length] = pvalue
                Corr_df.loc[start, length] = corr
                Datasize_df.loc[start, length] = len(x)
    print("======> P-values <=======")
    print(P_value_df)
    print("======> Spearman coor <=======")
    print(round(Corr_df,3))
    print("======> Datasize <=======")
    print(Datasize_df.astype(int))
    return P_value_df, Corr_df, Datasize_df

###################################
#### Read the TE data
###################################


mrna_05hr_1 = np.array([ d.replace(',','') if type(d)==str else d for d in media9.loc[:, 'mrna_05hr_1'].values ]).astype(float)
mrna_05hr_2 = np.array([ d.replace(',','') if type(d)==str else d for d in media9.loc[:, 'mrna_05hr_2'].values ]).astype(float)
mrna_24hr_1 = np.array([ d.replace(',','') if type(d)==str else d for d in media9.loc[:, 'mrna_24hr_1'].values ]).astype(float)
mrna_24hr_2 = np.array([ d.replace(',','') if type(d)==str else d for d in media9.loc[:, 'mrna_24hr_2'].values ]).astype(float)

fp_chx_05hr_1 = np.array([ d.replace(',','') if type(d)==str else d for d in media9.loc[:, 'fp_chx_05hr_1'].values ]).astype(float)
fp_chx_05hr_2 = np.array([ d.replace(',','') if type(d)==str else d for d in media9.loc[:, 'fp_chx_05hr_2'].values ]).astype(float)
fp_chx_24hr_1 = np.array([ d.replace(',','') if type(d)==str else d for d in media9.loc[:, 'fp_chx_24hr_1'].values ]).astype(float)
fp_chx_24hr_2 = np.array([ d.replace(',','') if type(d)==str else d for d in media9.loc[:, 'fp_chx_24hr_2'].values ]).astype(float)

condition1 = (mrna_05hr_1+mrna_05hr_2>10) & (mrna_24hr_1+mrna_24hr_2>10)
condition2 = (fp_chx_05hr_1+fp_chx_05hr_2>10) & (fp_chx_24hr_1+fp_chx_24hr_2>10)
condition = condition1 & condition2

mRNA_reads_05h = (mrna_05hr_1+mrna_05hr_2)[condition]
mRNA_reads_24h = (mrna_24hr_1+mrna_24hr_2)[condition]
RiboSeq_reads_05h = (fp_chx_05hr_1+fp_chx_05hr_2)[condition]
RiboSeq_reads_24h = (fp_chx_24hr_1+fp_chx_24hr_2)[condition]

names = media9['splice'].values[condition]

mRNA_reads_05h_dict = { k:v for k,v in zip(names, mRNA_reads_05h) }
mRNA_reads_24h_dict = { k:v for k,v in zip(names, mRNA_reads_24h) }
RiboSeq_reads_05h_dict = { k:v for k,v in zip(names, RiboSeq_reads_05h) }
RiboSeq_reads_24h_dict = { k:v for k,v in zip(names, RiboSeq_reads_24h) }

TE_05h = { k:RiboSeq_reads_05h_dict[k]/mRNA_reads_05h_dict[k] for k in set(mRNA_reads_05h_dict)&set(RiboSeq_reads_05h_dict) }
TE_24h = { k:RiboSeq_reads_24h_dict[k]/mRNA_reads_24h_dict[k] for k in set(mRNA_reads_24h_dict)&set(RiboSeq_reads_24h_dict) }

print(len(TE_05h))
print(len(TE_24h))


