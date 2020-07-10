
importCommon()
import SARS2

def read_UTR_icSHAPE():
    root = "/Share/home/zhangqf8/sunlei/data/SARS2/20200623-inVitro-SRAS-huh7-8UTR/8UTR-icSHAPE"
    icSHAPE = {}
    icSHAPE.update( General.load_shape(join(root, "MERS-HKU9/OUT/9.transSHAPE/virus.shape")) )
    icSHAPE.update( General.load_shape(join(root, "SARS2-C-NL63_HKU1/OUT/9.transSHAPE/virus.shape")) )
    icSHAPE.update( General.load_shape(join(root, "SARS2-T/OUT/9.transSHAPE/virus.shape")) )
    icSHAPE.update( General.load_shape(join(root, "SARS-HKU5/OUT/9.transSHAPE/virus.shape")) )
    return icSHAPE

utr_icSHAPE = read_UTR_icSHAPE()
utr_sequences = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/sequences/utr_sequences.fa')



#####################################
####  1. 首先通过调试预测出二级结构
#####################################


subseq = {}
subshape = {}
subbpprob = {}

subseq['SARS2-C'] = utr_sequences['SARS2-C'][-330:]
subshape['SARS2-C'] = utr_icSHAPE['SARS2-C'][-330:]
SARS2C_dot, subbpprob['SARS2-C'] = SARS2.maxexpect_predict(subseq['SARS2-C'], subshape['SARS2-C'])
cmd = Visual.Plot_RNAStructure_Shape(subseq['SARS2-C'], SARS2C_dot, subshape['SARS2-C'])

subseq['SARS'] = utr_sequences['SARS'][-330:]
subshape['SARS'] = utr_icSHAPE['SARS'][-330:]
SARS_dot, subbpprob['SARS'] = SARS2.maxexpect_predict(subseq['SARS'], subshape['SARS'])
cmd = Visual.Plot_RNAStructure_Shape(subseq['SARS'], SARS_dot, subshape['SARS'])


tmpseq1 = utr_sequences['MERS'][-290:-170]
tmpshape1 = utr_icSHAPE['MERS'][-290:-170]
tmpseq2 = utr_sequences['MERS'][-170:]
tmpshape2 = utr_icSHAPE['MERS'][-170:]
MERS_dot1, _ = SARS2.maxexpect_predict(tmpseq1, tmpshape1)
MERS_dot2, _ = SARS2.maxexpect_predict(tmpseq2, tmpshape2)
subseq['MERS'] = utr_sequences['MERS'][-290:]
subshape['MERS'] = utr_icSHAPE['MERS'][-290:]
MERS_dot = MERS_dot1 + MERS_dot2
cmd = Visual.Plot_RNAStructure_Shape(subseq['MERS'], MERS_dot, subshape['MERS'])
print(cmd)

tmpseq1 = utr_sequences['bHKU5'][-290:-170]
tmpshape1 = utr_icSHAPE['bHKU5'][-290:-170]
tmpseq2 = utr_sequences['bHKU5'][-170:]
tmpshape2 = utr_icSHAPE['bHKU5'][-170:]
bHKU5_dot1, _ = SARS2.maxexpect_predict(tmpseq1, tmpshape1)
bHKU5_dot2, _ = SARS2.maxexpect_predict(tmpseq2, tmpshape2)
subseq['bHKU5'] = utr_sequences['bHKU5'][-290:]
subshape['bHKU5'] = utr_icSHAPE['bHKU5'][-290:]
bHKU5_dot = bHKU5_dot1 + bHKU5_dot2
cmd = Visual.Plot_RNAStructure_Shape(subseq['bHKU5'], bHKU5_dot, subshape['bHKU5'])
print(cmd)


subseq['bHKU9'] = utr_sequences['bHKU9'][-265:]
subshape['bHKU9'] = utr_icSHAPE['bHKU9'][-265:]
bHKU9_dot, subbpprob['bHKU9'] = SARS2.maxexpect_predict(subseq['bHKU9'], subshape['bHKU9'])
cmd = Visual.Plot_RNAStructure_Shape(subseq['bHKU9'], bHKU9_dot, subshape['bHKU9'])
print(cmd)


subseq['hHKU1'] = utr_sequences['hHKU1'][-200:]
subshape['hHKU1'] = utr_icSHAPE['hHKU1'][-200:]
hHKU1_dot, subbpprob['hHKU1'] = SARS2.maxexpect_predict(subseq['hHKU1'], subshape['hHKU1'])
cmd = Visual.Plot_RNAStructure_Shape(subseq['hHKU1'], hHKU1_dot, subshape['hHKU1'])
print(cmd)


tmpseq1 = utr_sequences['hNL63'][-300:-250]
tmpshape1 = utr_icSHAPE['hNL63'][-300:-250]
tmpseq2 = utr_sequences['hNL63'][-250:]
tmpshape2 = utr_icSHAPE['hNL63'][-250:]
hNL63_dot1, _ = SARS2.maxexpect_predict(tmpseq1, tmpshape1)
hNL63_dot2, _ = SARS2.maxexpect_predict(tmpseq2, tmpshape2)
subseq['hNL63'] = utr_sequences['hNL63'][-300:]
subshape['hNL63'] = utr_icSHAPE['hNL63'][-300:]
hNL63_dot = hNL63_dot1 + hNL63_dot2
cmd = Visual.Plot_RNAStructure_Shape(subseq['hNL63'], hNL63_dot, subshape['hNL63'])
print(cmd)


#####################################
####  2. 收集数据
#####################################

subseq = {}
subshape = {}
subbpprob = {}

subseq['MERS'] = utr_sequences['MERS'][-290:]
subshape['MERS'] = utr_icSHAPE['MERS'][-290:]

subseq['SARS'] = utr_sequences['SARS'][-330:]
subshape['SARS'] = utr_icSHAPE['SARS'][-330:]

subseq['SARS2-C'] = utr_sequences['SARS2-C'][-330:]
subshape['SARS2-C'] = utr_icSHAPE['SARS2-C'][-330:]

subseq['SARS2-T'] = utr_sequences['SARS2-T'][-330:]
subshape['SARS2-T'] = utr_icSHAPE['SARS2-T'][-330:]

subseq['bHKU5'] = utr_sequences['bHKU5'][-290:]
subshape['bHKU5'] = utr_icSHAPE['bHKU5'][-290:]

subseq['bHKU9'] = utr_sequences['bHKU9'][-265:]
subshape['bHKU9'] = utr_icSHAPE['bHKU9'][-265:]

subseq['hHKU1'] = utr_sequences['hHKU1'][-200:]
subshape['hHKU1'] = utr_icSHAPE['hHKU1'][-200:]

subseq['hNL63'] = utr_sequences['hNL63'][-300:]
subshape['hNL63'] = utr_icSHAPE['hNL63'][-300:]

subbpprob = {}
for name in subseq:
    subbpprob[name] = Structure.partition(subseq[name], subshape[name], si=-0.4, sm=1.5, md=300)

structures = {
    'MERS':   '.....(((((....(((((.((...(((..(((.....)))...)))...)).)))))...)))))................(((((((((((...........))))))))))).........((((((((((..((((((...((.....(((((((..............((((((((...(((......))).))).))))).(((...))).......))))))).....))...))))))..)))))))))).((((((((......)))))).))........',
    'SARS':   '.......(((((((..(.(((((((((((..(((...((......))...))).)))))))))))))))))))...............(((((((((.((....))..)))))))))...((((..(((((.((((..((..(((.(((.....(((((((((.((.(((....))))).((.(((((((((....(((.((((.....))..))..)))...))))).))))...)).....((.((.....)).))..))))))))).....))).)))..))...)))).......)))))..))))....................',
    'SARS2-C':  '.......((((((..((.(((((((((((..(((...((......))...))).)))))))))))))))))))...............(((((((((.((....))..))))))))).((((((..((((.(((((..((..(((.(((.....(((((((((.((.(((....))))).((.(((((((((....((.((.((.........)).))))...))))).)))).((((........))))...)).....))))))))).....))).)))..))...))))).......))))..)))).)).................',
    'SARS2-T': '.......((((((..((.(((((((((((..(((...((......))...))).)))))))))))))))))))...............(((((((((.((....))..))))))))).((((((..((((.(((((..((..(((.(((.....(((((((((.((.(((....))))).((.(((((((((....((.((.((.........)).))))...))))).)))).((((........))))...)).....))))))))).....))).)))..))...))))).......))))..)))).)).................',
    'bHKU5':  '...(((((((...((((((.((..((((..(((.....)))..))))...)).)))))).)))))))................(((.(((((((...........))))))))))...........((((((((..((((......((.....(((((((..............((((((((...((((....)))).))).))))).(((...))).......))))))).....))......))))))))))))...(((..((((........))))..))).....',
    'bHKU9':  '.....(((((.(((.((((((((((..((((....))))..))))))).))).))).))))).................((((.((((((((((...........)))))))))).((((((.(((((((((((((...((.(((((...((((((((((((((................))).))..........)))))))))....))))).))...)))))))).....))))).))))))............))))....',
    'hHKU1':  '...((((((((((..(.....)..)))))))))).............((((..((((((((((.(((((((......))))....))).)))))))))).))))...........(((((((.(((.(((((.(((((..((.((((((((((((....))))))).)))))))))))))))))))))))))))......',
    'hNL63':  '......((((((((.(((((...(((....))))))))))))))))..........((((((((((((...)))))))....))))).........................(((((((((((..........)))))))))))(((((((.(((((((.......(((((((((((.....(((((((((((.............(((......)))((......)).....)))))))))))...))))).))))))))))).)).))))))).(((((......)))))........',
}

OUT = open(join(HOME, 'figs/icSHAPE-8UTR.ps1'), 'w')
for key in structures:
    auc = General.calc_AUC_v2(structures[key], subshape[key])
    title = f"{key} AUC={auc:.3}"
    cmd = Visual.Plot_RNAStructure_Shape(subseq[key], structures[key], subshape[key], title=title, mode='label', bpprob=subbpprob[key], bpprob_mode='both', bpwarning=False, )
    print(cmd, file=OUT)

OUT.close()

##########################
#### 寻找保守的碱基对
###########################

########  构建cm model

h = {}
for name in structures:
    stoFn = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_3UTR/{name}.sto'
    Covariation.dot2sto({'5UTR':[subseq[name], structures[name]]}, name, stoFn, mode='w')
    cmFn = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_3UTR/{name}.cm'
    Covariation.cmbuild(stoFn, cmFn, verbose=False, showCMD=True)
    h[name] = Covariation.cmcalibrate(cmFn, cpu=20, showCMD=True, use_LSF=True, LSF_parameters={})

########  Search from sequence

h = {}
for name in structures:
    seqdbFn = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/homoseq/Rfam_sequence.fa'
    outTXT = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_3UTR/cmsearch-{name}.txt'
    outSto = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_3UTR/cmsearch-{name}.sto'
    cmFn = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_3UTR/{name}.cm'
    h[name] = Covariation.cmsearch(cmFn, seqdbFn, outTXT, outSto, cpu=20, toponly=True, nohmm=True, nohmmonly=False, outputE=20, acceptE=1, verbose=True, showCMD=True, use_LSF=True, LSF_parameters={})

########  R-scape

for name in structures:
    outSto = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_3UTR/cmsearch-{name}.sto'
    Rscape_dir = '/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_3UTR/R-scape/'
    Covariation.R_scape(outSto, Rscape_dir, outname=name)

########  Calculate Covarying base pairs

bps = {}
for name in structures:
    clean_Sto = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_3UTR/R-scape/{name}.cons.sto'
    bps[name] = Covariation.calc_covBP_from_sto(clean_Sto, subseq[name], allpair=False, min_score=0.4)

########  Generate visualization command

def collect_covariate_bps(all_bps, dot):
    ct = Structure.dot2ct(dot)
    hl_regions = []
    for b1,b2,score in all_bps:
        if (b1,b2) in ct:
            if 0.4<=score<0.5:
                color = '#0877b5'#Colors.RGB['blue']
            elif 0.5<=score<0.7:
                color = '#c6900a'#Colors.RGB['orange']
            elif score>=0.7:
                color = '#09b509'#Colors.RGB['red']
            else:
                continue
            hl_regions.append([b1,b1,color])
            hl_regions.append([b2,b2,color])
        elif score>0.4:
            print(b1, b2, score)
    return hl_regions

cmd = {}
for name in structures:
    hl_regions = collect_covariate_bps(bps[name], structures[name])
    auc = General.calc_AUC_v2(structures[name], subshape[name])
    cmd[name] = Visual.Plot_RNAStructure_Shape(subseq[name], structures[name], 
        subshape[name], mode='label', bpprob=subbpprob[name],
        bpprob_cutofflist=[0.6, 0.8, 0.95], bpprob_mode='both', bpwarning=False,
        highlight_region=hl_regions, first_base_pos=1, period=10, peroid_color='#2196f3', 
        title=f"{name} AUC={auc:.3}")


OUT = open(join(HOME, 'figs/8UTR-3p.ps1'), 'w')
for name in cmd:
    print(cmd[name], file=OUT)

OUT.close()

