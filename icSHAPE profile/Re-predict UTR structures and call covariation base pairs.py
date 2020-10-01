import General, Structure, Visual, SARS2, Covariation

wuhan_id = 'NC_045512.2'

#shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus.shape')[wuhan_id]
shape = General.load_shape('/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/8.calcSHAPE/virus-w50.shape')[wuhan_id]

seq = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')[wuhan_id]
dot = General.load_dot('/Share2/home/zhangqf7/lipan/SARS2/2nd_structures/Rfam/structures.dot')

utr_5 = [1, 299]
utr_3 = [29536, 29870]
seq5 = seq[utr_5[0]-1:utr_5[1]]
shape5 = shape[utr_5[0]-1:utr_5[1]]
seq3 = seq[utr_3[0]-1:utr_3[1]]
shape3 = shape[utr_3[0]-1:utr_3[1]]

dot['RF03120_UTR5'][0]==seq5
dot['RF03125_UTR3'][0]==seq3

########  预测二级结构

prob3, pfs3 = Structure.partition(seq3, shape_list=shape3, si=-0.4, sm=1.5, md=300, return_pfs=True)
maxexpect_dot3 = Structure.maxExpect(input_pfs_file=pfs3, delete_pfs=True)[0]
clean_dot3 = SARS2.remove_lowporb_bp(prob3, maxexpect_dot3, minprob=0.6, minDist=150)
bp_prob3 = SARS2.collect_bpprob(prob3, clean_dot3)
auc3 = General.calc_AUC_v2(clean_dot3, shape3)
cmd3 = Visual.Plot_RNAStructure_Shape(seq3, clean_dot3, 
    shape3, mode='label', bpprob=bp_prob3,
    bpprob_cutofflist=[0.6, 0.8, 0.95], bpprob_mode='both',
    title=f"5UTR AUC={auc5:.3}")

prob5, pfs5 = Structure.partition(seq5, shape_list=shape5, si=-0.4, sm=1.5, md=300, return_pfs=True)
maxexpect_dot5 = Structure.maxExpect(input_pfs_file=pfs5, delete_pfs=True)[0]
clean_dot5 = SARS2.remove_lowporb_bp(prob5, maxexpect_dot5, minprob=0.6, minDist=150)
bp_prob5 = SARS2.collect_bpprob(prob5, clean_dot5)
auc5 = General.calc_AUC_v2(clean_dot5, shape5)
cmd5 = Visual.Plot_RNAStructure_Shape(seq5, clean_dot5, 
    shape5, mode='label', bpprob=bp_prob5,
    bpprob_cutofflist=[0.6, 0.8, 0.95], bpprob_mode='both',
    title=f"3UTR AUC={auc3:.3}")

########  构建cm model

utr5_sto = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/5UTR.sto'
Covariation.dot2sto({'5UTR':[seq5, clean_dot5]}, '5UTR', utr5_sto, mode='w')

utr3_sto = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/3UTR.sto'
Covariation.dot2sto({'3UTR':[seq3, clean_dot3]}, '3UTR', utr3_sto, mode='w')

utr5_cm = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/5UTR.cm'
Covariation.cmbuild(utr5_sto, utr5_cm, verbose=False, showCMD=True)

utr3_cm = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/3UTR.cm'
Covariation.cmbuild(utr3_sto, utr3_cm, verbose=False, showCMD=True)

h5 = Covariation.cmcalibrate(utr5_cm, cpu=20, showCMD=True, use_LSF=True, LSF_parameters={})
h3 = Covariation.cmcalibrate(utr3_cm, cpu=20, showCMD=True, use_LSF=True, LSF_parameters={})

########  Search from sequence

seqdbFn = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/homoseq/Rfam_sequence.fa'
outTXT5 = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/cmsearch-output-5.txt'
outSto5 = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/cmsearch-output-5.sto'
h5 = Covariation.cmsearch(utr5_cm, seqdbFn, outTXT5, outSto5, cpu=20, toponly=True, nohmm=True, nohmmonly=True, outputE=20, acceptE=1, verbose=True, showCMD=True, use_LSF=True, LSF_parameters={})

outTXT3 = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/cmsearch-output-3.txt'
outSto3 = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/cmsearch-output-3.sto'
h3 = Covariation.cmsearch(utr3_cm, seqdbFn, outTXT3, outSto3, cpu=20, toponly=True, nohmm=True, nohmmonly=True, outputE=20, acceptE=1, verbose=True, showCMD=True, use_LSF=True, LSF_parameters={})

########  R-scape filtering

outDir5 = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/5UTR_Rscape'
Covariation.R_scape(outSto5, outDir5, outname="5UTR")

outDir3 = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/3UTR_Rscape'
Covariation.R_scape(outSto3, outDir3, outname="3UTR")

########  Calculate Covarying base pairs

clean_Sto5 = "/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/5UTR_Rscape/5UTR.cons.sto"
all_bps5 = Covariation.calc_covBP_from_sto(clean_Sto5, seq5, allpair=True, min_score=0.2)

clean_Sto3 = "/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/3UTR_Rscape/3UTR.cons.sto"
all_bps3 = Covariation.calc_covBP_from_sto(clean_Sto3, seq3, allpair=True, min_score=0.2)

########  Generate visualization command

def collect_covariate_bps(all_bps, dot):
    ct = Structure.dot2ct(dot)
    hl_regions = []
    for b1,b2,score in all_bps:
        if (b1,b2) in ct:
            if 0.4<=score<0.5:
                color = Colors.RGB['blue']
            elif 0.5<=score<0.7:
                color = Colors.RGB['orange']
            elif score>=0.7:
                color = Colors.RGB['red']
            else:
                continue
            hl_regions.append([b1,b1,color])
            hl_regions.append([b2,b2,color])
        elif score>0.4:
            print(b1, b2, score)
    return hl_regions

hl_regions5 = collect_covariate_bps(all_bps5, clean_dot5)
hl_regions3 = collect_covariate_bps(all_bps3, clean_dot3)

cmd5 = Visual.Plot_RNAStructure_Shape(seq5, clean_dot5, 
    shape5, mode='label', bpprob=bp_prob5,
    bpprob_cutofflist=[0.6, 0.8, 0.95], bpprob_mode='both',
    highlight_region=hl_regions5, first_base_pos=1, period=10, peroid_color='#2196f3',
    title=f"5UTR AUC={auc5:.3}")
cmd3 = Visual.Plot_RNAStructure_Shape(seq3, clean_dot3, 
    shape3, mode='label', bpprob=bp_prob3,
    bpprob_cutofflist=[0.6, 0.8, 0.95], bpprob_mode='both',
    highlight_region=hl_regions3, first_base_pos=29536, period=10, peroid_color='#2196f3',
    title=f"3UTR AUC={auc3:.3}")



