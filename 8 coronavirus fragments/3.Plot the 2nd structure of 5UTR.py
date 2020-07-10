
importCommon()
import SARS2
import Covariation

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

def refine_structure(seq1, dot1, seq2, shape2):
    """
    Use seq1 and dot1 to refine seq2
    """
    alignments = Structure.multi_alignment([seq1, seq2])
    aligned_dot = Structure.dot_to_alignDot(dot1, alignments[0])
    stem_loops = SARS2.find_stemloop(aligned_dot)
    
    dot2 = "."*len(seq2)
    for sl in stem_loops:
        sub_aligned_seq = alignments[1][sl[0]-1:sl[1]]
        sub_seq = sub_aligned_seq.replace('-', '')
        start = seq2.find(sub_seq)
        end = start+len(sub_seq)
        sub_dot = Structure.predict_structure(sub_seq, shape2[start:end])
        dot2 = dot2[:start]+sub_dot+dot2[end:]
    
    return dot2


subseq = {}
subshape = {}
subbpprob = {}

subseq['MERS'] = utr_sequences['MERS'][:460]
subshape['MERS'] = utr_icSHAPE['MERS'][:460]

subseq['SARS'] = utr_sequences['SARS'][:395]
subshape['SARS'] = utr_icSHAPE['SARS'][:395]

subseq['SARS2-C'] = utr_sequences['SARS2-C'][:395]
subshape['SARS2-C'] = utr_icSHAPE['SARS2-C'][:395]

subseq['SARS2-T'] = utr_sequences['SARS2-T'][:395]
subshape['SARS2-T'] = utr_icSHAPE['SARS2-T'][:395]

subseq['bHKU5'] = utr_sequences['bHKU5'][:445]
subshape['bHKU5'] = utr_icSHAPE['bHKU5'][:445]

subseq['bHKU9'] = utr_sequences['bHKU9'][:330]
subshape['bHKU9'] = utr_icSHAPE['bHKU9'][:330]

subseq['hHKU1'] = utr_sequences['hHKU1'][:330]
subshape['hHKU1'] = utr_icSHAPE['hHKU1'][:330]

subseq['hNL63'] = utr_sequences['hNL63'][:330]
subshape['hNL63'] = utr_icSHAPE['hNL63'][:330]

MERS_dot, subbpprob['MERS'] = SARS2.maxexpect_predict(subseq['MERS'], subshape['MERS'])
SARS_dot, subbpprob['SARS'] = SARS2.maxexpect_predict(subseq['SARS'], subshape['SARS'])
SARS2C_dot, subbpprob['SARS2-C'] = SARS2.maxexpect_predict(subseq['SARS2-C'], subshape['SARS2-C'])
SARS2T_dot, subbpprob['SARS2-T'] = SARS2.maxexpect_predict(subseq['SARS2-T'], subshape['SARS2-T'])
subbpprob['bHKU5'] = Structure.partition(subseq['bHKU5'], subshape['bHKU5'], si=-0.4, sm=1.5, md=300)
MERS_str = ".....(((((((((((...))))).))))))......(((((.....))))).(((.......)))............((((.(.((((.(((((((.(((.((((.((((((((((.....((((......))))..)))))).)))))))))))))))))).))))).)))).....((.....)).......(((((((...((((((..((.(((..(((((((((((((((..(((.(((......)))))).)))))....)))).(((((((.(((......))))))))))(((((((.......))))))))))))).)))))))))))....))))))).......(((((.((.((.....(((.(((((...))))).))))).)).)))))......((((((((((..((((((...((((....))))))))))))))))))))."
bHKU5_dot = refine_structure(subseq['MERS'], MERS_str, subseq['bHKU5'], subshape['bHKU5'])
bHKU9_dot, subbpprob['bHKU9'] = SARS2.maxexpect_predict(subseq['bHKU9'], subshape['bHKU9'])
hHKU1_dot, subbpprob['hHKU1'] = SARS2.maxexpect_predict(subseq['hHKU1'], subshape['hHKU1'])
hNL63_dot, subbpprob['hNL63'] = SARS2.maxexpect_predict(subseq['hNL63'], subshape['hNL63'])

structures = {
    'MERS':   '.....(((((((((((...))))).))))))......(((((.....))))).(((.......)))............((((.(.((((.(((((((.(((.((((.((((((((((.....((((......))))..)))))).)))))))))))))))))).))))).)))).....................(((((((...((((((..((.(((..(((((((((((((((..(((.(((......)))))).)))))....)))).(((((((.(((......))))))))))(((((((.......))))))))))))).)))))))))))....))))))).......(((((.((.((.....(((.(((((...))))).))))).)).)))))......((((((((((..((((((...((((....)))))))))))))))))))).',
    'SARS':   '......(((((((.(((....))))))))))..........(((((.....))))).((((.......))))........((((((((.((.((((.(((.....))).)))))).))))))))........................((((((((((((.(((((...(((.(((.(((((((..((((((.(((((......)))))..))))))......)))(((((((.((......)))))))))(((....))))))).)))))).))))))))))...))))))).......((((((....((.....(((((.....)))))..)))))))).....(((((.((((((((.((((.....)))).)))...))))).)))))..',
    'SARS2-C':  '......(((((.(((((....)))))..)))))...........(((((.....))))).((((.......))))........((((((((.((.((((.(((.....))).)))))).))))))))......................(((((((((((..(((((...(((.(((((((((((..((((((.(((((......)))))..))))))......)))(((((((.((......)))))))))(((....)))))))))))))).))))).))))...))))))).......((((((.(((...))).((((((...))))))....)))))).....(((((.(((((((((((((.....)))).))))..))))).))))).',
    'SARS2-T': '......(((((.(((((....)))))..)))))...........(((((.....))))).((((.......))))........((((((((.((.((((.(((.....))).)))))).))))))))......................(((((((((((..(((((...(((.(((((((((((..((((((.(((((......)))))..))))))......)))(((((((.((......)))))))))(((....)))))))))))))).))))).))))...))))))).......((((((.(((...))).((((((...))))))....)))))).....(((((.(((((((((((((.....)))).))))..))))).))))).',
    'bHKU5':  '......((((.(((((...)))))..)))).....(((((.....))))).(((((...)))))............((((...((((.(((((((.((....((..((((..(((((...)))))))))..))...))))))))).))))..)))).....................((((((((..(((..((.((.(((..(((((((((((((((..(((.(((......)))))).)))))....).)))(((((((.(((......))))))))))(((((((.......))))))))))))).))))))).)))...)))))))).......(((((...((((((....(((((((...))))))).)))))).)))))........(((((.(((.((((((...((((....))))))))))))))))))......',
    'bHKU9':  '.......(((.(((((.....))))).)))...............(((((.....)))))(((....))).............((((...((((.(((((((((......)))))))))))))...))))..............((((((...(((((.((((..(((((..((.(((((((((......))))).)))).)).((((((((......))))))))....)).))).)))).))))).....))))))..............(((((((....)))))))..............((((((((.......))))))))...',
    'hHKU1':  '.....((((((((.((((.......)))).))).))))).(((((......)))))......................((((((((.........(((((((....))))))).........))))))))..................(((((((....(((.((....((((((((....((((((.((........)).)))..)))....))))))))...)))))...(((((((...(((((((((...(((.((.(((((........))))).)).))))))))..))))..)))))))))))))).................',
    'hNL63':  '....(((((.(((((((((....)))))))))..)))))...(((((.....)))))..........((..(((((((((..((((.((((((((.....)))))))).))))........)))))))))..)).(((((.((((((((.(((.((((((.......(((((....((((((((((......)))))))))))))))........(((((((((......)))))))))..((....((((((......))))))....))))))).)))))))))).))..))))).....((((((((((....)))))))..)))..'
}

OUT = open(join(HOME, 'figs/icSHAPE-8UTR.ps1'), 'w')
for key in structures:
    auc = General.calc_AUC_v2(structures[key], subshape[key])
    title = f"{key} AUC={auc:.3}"
    cmd = Visual.Plot_RNAStructure_Shape(subseq[key], structures[key], subshape[key], title=title, mode='label', bpprob=subbpprob[key], bpprob_mode='both' )
    print(cmd, file=OUT)

OUT.close()


###########################
#### 寻找保守的碱基对
###########################

########  构建cm model

h = {}
for name in structures:
    stoFn = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_5UTR/{name}.sto'
    Covariation.dot2sto({'5UTR':[subseq[name], structures[name]]}, name, stoFn, mode='w')
    cmFn = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_5UTR/{name}.cm'
    Covariation.cmbuild(stoFn, cmFn, verbose=False, showCMD=True)
    h[name] = Covariation.cmcalibrate(cmFn, cpu=20, showCMD=True, use_LSF=True, LSF_parameters={})

########  Search from sequence

h = {}
for name in structures:
    seqdbFn = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/homoseq/Rfam_sequence.fa'
    outTXT = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_5UTR/cmsearch-{name}.txt'
    outSto = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_5UTR/cmsearch-{name}.sto'
    cmFn = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_5UTR/{name}.cm'
    h[name] = Covariation.cmsearch(cmFn, seqdbFn, outTXT, outSto, cpu=20, toponly=True, nohmm=True, nohmmonly=False, outputE=20, acceptE=1, verbose=True, showCMD=True, use_LSF=True, LSF_parameters={})

# /Share/home/zhangqf8/usr/infernal/bin/cmsearch --notextw --cpu 40 --toponly --nohmm -E 20 --incE 10 -o /Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_5UTR/cmsearch-bHKU9.txt -A /Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_5UTR/cmsearch-bHKU9.sto /Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_5UTR/bHKU9.cm /Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/homoseq/Rfam_sequence.fa

########  R-scape

for name in structures:
    outSto = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_5UTR/cmsearch-{name}.sto'
    Rscape_dir = '/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_5UTR/R-scape/'
    Covariation.R_scape(outSto, Rscape_dir, outname=name)

########  Calculate Covarying base pairs

bps = {}
for name in structures:
    clean_Sto = f'/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/call_covariation_5UTR/R-scape/{name}.cons.sto'
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
        subshape[name], mode='label', bpprob=subbpprob[name], bpwarning=False,
        bpprob_cutofflist=[0.6, 0.8, 0.95], bpprob_mode='both',
        highlight_region=hl_regions, first_base_pos=1, period=10, peroid_color='#2196f3',
        title=f"{name} AUC={auc:.3}")

OUT = open(join(HOME, 'figs/8UTR-5p.ps1'), 'w')
for name in cmd:
    print(cmd[name], file=OUT)

OUT.close()



