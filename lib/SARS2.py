
from Bio import AlignIO
from Bio.Alphabet import generic_rna
import Structure, General

wuhan_id = 'NC_045512.2'
virus_seq = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')[wuhan_id]
virus_shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-w50.shape')[wuhan_id]
virus_dot = General.load_dot('/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/SARS2.dot')[wuhan_id][1]

def remove_lowporb_bp(prob_list, dot, minprob=0.6, minDist=100):
    """
    prob_list:          [(1,19,0.99211), ...]
    dot:                Dot-Bracket
    minprob/minDist     If a base pair longer than minDist with probability less than minprob,
                        it will be removed
    """
    ctlist = Structure.dot2ct(dot)
    nCtList = []
    i,j = 0,0
    while i<len(ctlist) and j<len(prob_list):
        if ctlist[i][1]-ctlist[i][0]<minDist:
            nCtList.append( ctlist[i] )
        else:
            while j<len(prob_list) and ctlist[i]>prob_list[j][:2]:
                j += 1
            if j<len(prob_list):
                if ctlist[i]!=prob_list[j][:2]:
                    print( Colors.f(f"Warning: base pair {ctlist[i][0]}-{ctlist[i][1]} is not a member of prob_list", fc='yellow') )
                else:
                    if prob_list[j][2]>=minprob:
                        nCtList.append( ctlist[i] )
        i += 1
    nDot = Structure.ct2dot(nCtList,len(dot))
    return nDot

def find_stemloop(dot):
    stemloop = []
    stack = []
    for i,label in enumerate(dot):
        if label=='(':
            stack.append( i )
        elif label==')':
            if len(stack)==1:
                lb = stack[0]
                stemloop.append((lb+1, i+1))
            stack.pop()
    assert len(stack)==0, "stack size should be empty!"
    return stemloop


def collect_bpprob(prob_list, dot):
    """
    prob_list               -- prob_list:          [(1,19,0.99211), ...]
    dot                     -- ....((((....))))..
    
    Get subset of pairing probability from prob_list in dot
    """
    bp_prob = []
    ct_list = Structure.dot2ct(dot)
    i = 0
    j = 0
    while i<len(prob_list) and j<len(ct_list):
        while i<len(prob_list) and prob_list[i][:2]<ct_list[j]:
            i += 1
        if i<len(prob_list) and prob_list[i][:2]==ct_list[j]:
            bp_prob.append( prob_list[i] )
            i += 1
            j += 1
        elif i<len(prob_list):
            j += 1
    return bp_prob

def dot2sto(dot, modelname, outfile, mode='w'):
    """
    dot             -- { seqname:[aligned_seq, aligned_dot], ... }
    modelname       -- CM name
    outfile         -- Write the model to file
    mode            -- Cover or append
    """
    OUT = open(outfile, mode)
    print("# STOCKHOLM 1.0\n", file=OUT)
    print(f"#=GF ID   {modelname}\n", file=OUT)
    common_dot = ""
    for seqName in dot:
        align_seq, align_dot = dot[seqName]
        if not common_dot:
            common_dot = align_dot
        else:
            assert align_dot==common_dot, "The aligned dot should be same"
        print("%-20s%s" % (seqName, align_seq), file=OUT)
    print("%-20s%s" % ("#=GC SS_cons", common_dot), file=OUT)
    print("\n//", file=OUT)
    OUT.close()

def load_stockholm(stoFn):
    stoObj = AlignIO.read(stoFn, "stockholm")
    alignObjs = list(iter(stoObj))
    id2seq = { alignObj.id:str(alignObj.seq) for alignObj in alignObjs }
    structure = stoObj.column_annotations.get("secondary_structure", "")
    refA = stoObj.column_annotations.get("reference_annotation", "")
    return id2seq, structure, refA

def read_Rscape(Rscape_sorted_cov_fn):
    cov_pairs = []
    for line in open(Rscape_sorted_cov_fn):
        if line[0]=='*':
            data = line.strip().split()
            left = int(data[1])
            right = int(data[2])
            cov_pairs.append((left, right))
    return sorted(cov_pairs, key=lambda x: x[0])

def get_alignedPos2cleanPos_dict(refA):
    refA = refA.upper().replace('~','-').replace('.','-').replace(':','-').replace(',','-')
    alignedPos2cleanPos = {}
    i, j = 0, 0
    while i<len(refA):
        if refA[i]=='-':
            alignedPos2cleanPos[i+1] = None
        else:
            alignedPos2cleanPos[i+1] = j+1
            j += 1
        i += 1
    return alignedPos2cleanPos

def collect_bpbases(Rscape_inSto, Rscape_usedSto, raw_pos1, raw_pos2):
    import collections
    id2seq, structure, refA = load_stockholm(Rscape_inSto)
    used_keys = load_stockholm(Rscape_usedSto)[0].keys()
    used_keys = list(used_keys)
    for i in range(len(used_keys)):
        data = used_keys[i].split('_')
        used_keys[i] = "_".join(data[:-1]) + "/" + data[-1]
    rscape_list = read_Rscape(rscape_file)
    refA = refA.upper().replace('~','-').replace('.','-').replace(':','-').replace(',','-')
    alignedPos2cleanPos = get_alignedPos2cleanPos_dict(refA)
    cleanPos2alignedPos = { v:k for k,v in alignedPos2cleanPos.items() }
    
    aligned_pos1 = cleanPos2alignedPos[raw_pos1]
    aligned_pos2 = cleanPos2alignedPos[raw_pos2]
    
    bp_list = []
    for key in used_keys:
        b1 = id2seq[key][aligned_pos1-1]
        b2 = id2seq[key][aligned_pos2-1]
        bp_list.append( b1+b2 )
    count = collections.Counter(bp_list)
    return dict(count)

def print_covariation_sites(stoFn, Rscape_file, Rscape_usedSto):
    """
    Print Covariation Sites
    stoFn                       -- Stockholm file (Rscape input sto file)
    Rscape_file                 -- Rscape covariation file (.sorted.cov)
    Rscape_usedSto              -- Rscape output sto file (.cons.sto)
    """
    id2seq, structure, refA = load_stockholm(stoFn)
    rscape_list = read_Rscape(Rscape_file)
    
    refA = refA.upper().replace('~','-').replace('.','-').replace(':','-').replace(',','-')
    alignedPos2cleanPos = get_alignedPos2cleanPos_dict(refA)
    
    renumed_rscape_list = []
    for b1,b2 in rscape_list:
        renumed_rscape_list.append( (alignedPos2cleanPos[b1], alignedPos2cleanPos[b2]) )
    
    cleanSeq = refA.replace('-', '')
    for b1,b2 in renumed_rscape_list:
        ls = max(b1-1-5, 0)
        le = min(b1+5, len(cleanSeq))
        rs = max(b2-1-5, 0)
        re = min(b2+5, len(cleanSeq))
        left_seq = cleanSeq[ls:b1-1]+Colors.f(cleanSeq[b1-1], fc='red')+cleanSeq[b1:le]
        right_seq = cleanSeq[rs:b2-1]+Colors.f(cleanSeq[b2-1], fc='red')+cleanSeq[b2:re]
        count = collect_bpbases(stoFn, Rscape_usedSto, b1, b2)
        print(f"=======> {b1} {b2} <========", count)
        print( left_seq, f"{ls}-{le}" )
        print( right_seq, f"{rs}-{re}" )

def maxexpect_predict(subseq, subshape, si=-0.4, sm=1.5, md=300, remove_lowbp=True):
    prob, pfs = Structure.partition(subseq, shape_list=subshape, si=si, sm=sm, md=md, verbose=False, return_pfs=True)
    dot = Structure.maxExpect(input_pfs_file=pfs, delete_pfs=True)[0]
    # Remove base pairs long than 150nt with pairing probability less than 60%
    if remove_lowbp:
        dot = remove_lowporb_bp(prob, dot, minprob=0.6, minDist=150)
    bpprob = collect_bpprob(prob, dot)
    return dot, bpprob

def call_covariation_split(fullseq, fulldot, seqdbFn, workdir_root, nohmm=True):
    import Covariation, os
    
    all_covary_bps = []
    sl_list = find_stemloop(fulldot)
    i = 0
    while i<len(sl_list):
        start, end = sl_list[i]
        while i<len(sl_list)-1 and end-start<50:
            i += 1
            end = sl_list[i][1]
        i += 1
        print(f"====>>> {start}-{end} {fulldot[start-1:end]}")
        model_name = f"{start}-{end}"
        workdir = os.path.join(workdir_root, model_name)
        covary_bps = Covariation.call_covariation(fullseq[start-1:end], fulldot[start-1:end], 
            model_name, seqdbFn, workdir=workdir,
            nohmm=nohmm, cmsearchE=1, cpu=20, use_LSF=True, 
            LSF_parameters={}, progress=True, clean=False)
        for bp in covary_bps:
            all_covary_bps.append((bp[0]+start-1, bp[1]+start-1))
    return all_covary_bps


def read_ORF():
    # /Share2/home/zhangqf7/Jbrowser/SARS2/SARS2.gff
    ORF = {
        '5\'UTR': [1, 265],
        'ORF1ab': [266, 21555],
        'S': [21563, 25384],
        'ORF3a': [25393, 26220],
        'E': [26245, 26472],
        'M': [26523, 27191],
        'ORF6': [27202, 27387],
        'ORF7a': [27394, 27759],
        'ORF7b': [27756, 27887],
        'ORF8': [27894, 28259],
        'N': [28274, 29533],
        'ORF10': [29558, 29674],
        '3\'UTR': [29675, 29903]
    }
    return ORF

def annotate_region(start, end):
    region_list = []
    ORF = read_ORF()
    for name in ORF:
        if ORF[name][0]<end and start<ORF[name][1]:
            region_list.append(name)
    return " ".join(region_list)

#stoFn = "/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/3UTR_search.sto"
#Rscape_file = "/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/3UTR_Rscape/3UTR_search_1.sorted.cov"
#Rscape_usedSto = "/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/3UTR_Rscape/3UTR_search_1.cons.sto"
#print_covariation_sites(stoFn, Rscape_file, Rscape_usedSto)


#Rscape_inSto = "/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/3UTR_search.sto"
#Rscape_usedSto = "/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/3UTR_Rscape/3UTR_search_1.cons.sto"
#collect_bpbases(Rscape_inSto, Rscape_usedSto, 193, 233)
#collect_bpbases(Rscape_inSto, Rscape_usedSto, 193+1, 233-1)
#collect_bpbases(Rscape_inSto, Rscape_usedSto, 193+2, 233-2)
#collect_bpbases(Rscape_inSto, Rscape_usedSto, 193+3, 233-3)






#id2seq, structure = load_stockholm("/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/5UTR_search.sto")
#id2seq_used, structure_used = load_stockholm("/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/UTR_fitting/5UTR_Rscape/5UTR_search_1.cons.sto")

#p1 = 327
#p2 = 357
#for key in id2seq_used:
#    key = key.split('_')
#    key = "_".join(key[:-1])+"/"+key[-1]
#    s1 = id2seq[key][p1-1]
#    s2 = id2seq[key][p2-1]
#    if key.startswith('NC_045512.2'):
#        print("==============>", key)
#    if s1!='-':
#        print(s1, s2, id2seq[key][p1-1-12:p1-1+12], id2seq[key][p2-1-12:p2-1+12], structure[p1-1-12:p1-1+12], structure[p2-1-12:p2-1+12], key)




#dot  = "::::::::::::<<<<<<<--<-<<<<<<<<<<<--<<<---<<______>>--->>>.->>>>>>>>>>>>>>>>>>>--------{{{{{{,.<<<<<<<<<-<<____>>-->>>>>>>>>,[[[[[[--[[[[[[[[[[--.[[--[[[-[.[[-----[[[[[[[[[,<<-<<<____>..>>>>,((,~~~~~~~~~~~~,,,)),,,,,]]]]]]]]]-----]].]-]]]--]]---]]]]]------]]]]]--]]]]-]]<<____>>}}}}}}:::"
#seq1 = "UCAUGCAGACCACaCaaGgCAGAugGgCuauauaAACGuUUUCGCUUUUCCGUUUaCG.AuauauaGcCcaCcCuuGuGCAGAAUGAauuCuCG.uaaCuaCauAgCACAAGcAGauGuaGuuaACuuuaaUCuCaCauaGCaAU.CuUUaauCa.guGUGUAaCauuaGGGAGGACucGAAAg..aGCCACCA*[50]**[17]*UAUGGAAGAGCCCuaauGuGUAAAac.uAauuUUaGUAGuGCuaUCCCCAuGuGaUUuuaaUaGCuUCUUaGGaGaauGAC"
#seq2 = "UCACUCAAAGUA-ACAAGAUCGC--GGCAAUCGUUUGUGUUUGGCAACCCCAUCUCACcAUCGCUUGUCCACUCUUGCACAGAAUGGAAUCAUGuUGUAAUUACAGUGCAAUAAGGUAAUUAUAACCCAUUUAAUUGAUAGCUAUgCUUUAUUAAaGUGUGUAGCUGUAGAGAGAAUGUUAAAGacUGUCACCU*[13]**[18]*---GGAAGAGCUCUACAGUGUGAAAUgUAAAUAAAAAAUAGCUAUUAU--UCAAUUAGAUUAGGCUAAUUAGAUGAUUUGC"
#ct = Structure.dot2ct(dot)

#new_ct = []

#valid = 0
#invalid = 0
#for b1,b2 in ct:
#    if seq2[b1-1].upper()+seq2[b2-1].upper() in ('AU','UA','GC','CG','GU','UG'):
#        valid += 1
#        new_ct.append( (b1,b2) )
#    else:
#        print( seq2[b1-1].upper()+seq2[b2-1].upper() )
#        invalid += 1

#print(valid, invalid)


#new_dot = Structure.ct2dot(new_ct, len(seq2))
#new_dot2 = ""
#new_seq2 = ""
#for s,d in zip(seq2, new_dot):
#    if s.upper() in ('A','U','C','G','T'):
#        new_seq2 += s.upper()
#        new_dot2 += d

