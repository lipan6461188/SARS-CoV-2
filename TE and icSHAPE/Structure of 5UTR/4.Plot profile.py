
importCommon()

def readRTBD(gTabFn, size):
    RT = {}
    BD = {}
    for key in size:
        RT[key] = [[0,0] for _ in range(size[key])]
        BD[key] = [[0,0] for _ in range(size[key])]
    for line in open(gTabFn):
        if line[0]!='@':
            data = line.strip().split()
            if data[1]=='-':
                continue
            key = data[0]
            pos = int(data[2])
            RT[key][pos-1] = (int(data[3]), int(data[5]))
            BD[key][pos-1] = (int(data[4]), int(data[6]))
    return RT, BD


splice_gTab = "/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/13.map_virus_splice/only_cover_splice/shape.gTab"
full_gTab = "/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/8.calcSHAPE/virus.gTab"
splice_chrNameLength = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/splice_fasta/chrNameLength.txt'
full_chrNameLength = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/index/SARS2_complement.size'
full_size = { line.strip().split()[0]:int(line.strip().split()[1]) for line in open(full_chrNameLength) }
splice_size = { line.strip().split()[0]:int(line.strip().split()[1]) for line in open(splice_chrNameLength) }

splice_RT, splice_BD = readRTBD(splice_gTab, splice_size)
full_RT, full_BD = readRTBD(full_gTab, full_size)

keys1 = ['63_21549', '65_21551', '65_25380', '69_26236', '64_26467', '69_27040', '66_27384', '65_27883', '64_28254', '68_28262']
keys2 = ['NC_045512.2']

DMSO_RT = { 'full':[d[1] for d in full_RT['NC_045512.2']] }
DMSO_BD = { 'full':[d[1] for d in full_BD['NC_045512.2']] }
NAI_RT =  { 'full':[d[0] for d in full_RT['NC_045512.2']] }
NAI_BD =  { 'full':[d[0] for d in full_BD['NC_045512.2']] }
for key in keys1:
    NAI_RT[key] =  [d[0] for d in splice_RT[key]]
    DMSO_RT[key] = [d[1] for d in splice_RT[key]]
    NAI_BD[key] =  [d[0] for d in splice_BD[key]]
    DMSO_BD[key] = [d[1] for d in splice_BD[key]]


def plot_profile(dataset, savefn):
    plt.figure(figsize=(12, 15))
    for ii,key in enumerate(dataset):
        data = dataset[key][:65]
        data = [ 0 if d=='NULL' else float(d) for d in data ]
        plt.subplot(11, 1, ii+1)
        plt.bar(range(1, 66), data)
        plt.title(key)
    plt.tight_layout()
    plt.savefig(savefn)
    plt.close()

plot_profile(NAI_RT, join(HOME, 'figs/NAI_rt.pdf'))
plot_profile(DMSO_RT, join(HOME, 'figs/DMSO_rt.pdf'))
plot_profile(NAI_BD, join(HOME, 'figs/NAI_bd.pdf'))
plot_profile(DMSO_BD, join(HOME, 'figs/DMSO_bd.pdf'))


########################
### Shape profile
########################


root = "/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing"
fullshape = General.load_shape(join(root, "8.calcSHAPE/virus-w50.shape"))
#spliceshape = General.load_shape(join(root, "13.map_virus_splice/only_cover_splice/final.shape"))
spliceshape = General.load_shape(join(root, "14.map_virus_splice-all/only_cover_splice/final.shape"))

shape_dict = {'full': fullshape['NC_045512.2']}
#keys1 = ['63_21549', '65_21551', '65_25380', '69_26236', '64_26467', '69_27040', '66_27384', '65_27883', '64_28254', '68_28262']
keys1 = ['68_18050', '64_28254']
for key in keys1:
    shape_dict[key] = spliceshape[key]


def plot_profile(dataset, savefn):
    plt.figure(figsize=(12, 15))
    for ii,key in enumerate(dataset):
        if key!='full':
            end = int(key.split('_')[0])
        else:
            end = 10000
        raw_data = dataset[key][:68]
        data = [ 1 if d=='NULL' or i>=end else float(d) for i,d in enumerate(raw_data) ]
        colors = [ Colors.RGB['gray'] if d=='NULL' or i>=end else Colors.RGB['blue'] for i,d in enumerate(raw_data) ]
        plt.subplot(11, 1, ii+1)
        plt.bar(range(1, 69), data, color=colors)
        plt.title(key)
    plt.tight_layout()
    plt.savefig(savefn)
    plt.close()

plot_profile(shape_dict, join(HOME, 'figs/icSHAPE_profile.pdf'))





