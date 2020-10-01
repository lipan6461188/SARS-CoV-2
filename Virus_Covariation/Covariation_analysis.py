 
# cm <==> seqDB
# cmsearch => 去重等 => cmalign => cmbuild => cmcalibrate => cmsearch ...


# input_sto = '/Share2/home/zhangqf7/tmp/covariation/SASR2_validation/148-296/input.sto'
# input_seq = 'UGUCGUUGACAGGACACGAGUAACUCGUCUAUCUUCUGCAGGCUGCUUACGGUUUCGUCCGUGUUGCAGCCGAUCAUCAGCACAUCUAGGUUUCGUCCGGGUGUGACCGAAAGGUAAGAUGGAGAGCCUUGUCCCUGGUUUCAACGAGA'
# seqDBfn = '/Share2/home/zhangqf7/tmp/covariation/GenomicFastaResults.uniq.fa'
# workdir = '/Share2/home/zhangqf7/tmp/covariation/SASR2_validation/148-296'
# iterate_alignstr(input_sto, input_seq, seqDBfn, workdir, iteration=3)

def iterate_alignstr(input_sto, input_seq, seqDBfn, workdir, iteration=3, blastFilter=None):
    """
    blastFilter             -- If set, will use blast to filter the cmsearch target
                               { 'full_seq': "ACTG...", 
                                'blastdb':"", 
                                'flanking': 20, 
                                'blast_identity': 0.4,
                                'ignore_seq': ['input'] }
    """
    
    if blastFilter:
        import Alignment
        input_start = blastFilter['full_seq'].find(input_seq)
        assert input_start != -1, f"input_seq not in full_seq"
        input_end = input_start + len(input_seq)
        start = max(input_start-blastFilter.get('flanking', 20), 0)
        end = min(input_end+blastFilter.get('flanking', 20), len(blastFilter['full_seq']))
        hits = Alignment.blast_sequence_V2(blastFilter['full_seq'][start:end], blastFilter['blastdb'], verbose=True, task='blastn', maxhit=5000, perc_identity=blastFilter['blast_identity'] )
        target_valid_range = {}
        for hit in hits:
            if hit.hit_acc not in target_valid_range:
                target_valid_range[hit.hit_acc] = (hit.hit_from, hit.hit_to)
        if len(target_valid_range)<5:
            print(Colors.f("Too little blast target"))
            return False
        else:
            print(f"Blast: {len(target_valid_range)} results found")
    
    for it_num in range(iteration):
        print(f"===========>>> Iteration {it_num}  <<<===========")
        
        ### Define files 
        if it_num==0:
            input_sto = input_sto
        else:
            input_sto = join(workdir, f'cmalign_output_{it_num-1}.sto')
        input_cm = join(workdir, f'cmsearch_input_{it_num}.cm')
        cmsearch_output_txt = join(workdir, f'cmsearch_output_{it_num}.txt')
        cmsearch_output_sto = join(workdir, f'cmsearch_output_{it_num}.sto')
        cmalign_input_fa = join(workdir, f'cmalign_input_{it_num}.fa')
        cmalign_output_sto = join(workdir, f'cmalign_output_{it_num}.sto')
        
        ### cmsearch
        Covariation.cmbuild(input_sto, input_cm)
        Covariation.cmcalibrate(input_cm, use_LSF=True, LSF_parameters={'cpu': 2}).wait()
        tmp_fa,tmp_annot = General.load_fasta(seqDBfn, load_annotation=True)
        tmp_fa['input'] = input_seq
        General.write_fasta(tmp_fa, join(workdir, 'tmp_ref.fa'), tmp_annot)
        del tmp_fa
        h = Covariation.cmsearch(input_cm, join(workdir, 'tmp_ref.fa'), \
            cmsearch_output_txt, cmsearch_output_sto, cpu=20, toponly=True, nohmm=False, nohmmonly=False, \
            outputE=20, acceptE=1, cut_ga=False, rfam=False, glocal=False, verbose=True, showCMD=True, use_LSF=True, 
            LSF_parameters={'cpu': 5})
        h.wait()
        os.remove(join(workdir, 'tmp_ref.fa'))
        
        ### cmalign
        id2seq_dict, refStr, refSeq = General.load_stockholm(cmsearch_output_sto)[0]
        refID = [ key for key in id2seq_dict if key.startswith('input') ][0]
        refSeq = id2seq_dict[refID]
        
        ### 过滤不在blast结果中的区域
        #gb:LR821888|Organism:Severe/482-557
        if blastFilter:
            for key in list(id2seq_dict.keys()):
                if '/' in key:
                    true_id, start_end = key.rsplit('/', 1)
                    start, end = start_end.split('-')
                    start = int(start)
                    end = int(end)
                    if true_id in target_valid_range:
                        valid_start, valid_end = target_valid_range[true_id]
                        if (end<valid_start or start>valid_end):
                            del id2seq_dict[key]
                    elif true_id in blastFilter.get('ignore_seq',['input']):
                        print("Ignored:",true_id)
                    else:
                        print("Killed:",Colors.f(true_id, fc='red'))
                        del id2seq_dict[key]
                else:
                    print("Killed:",Colors.f(true_id, fc='green'))
                    del id2seq_dict[key]
        
        ### 去除相似的的序列
        uniq_seq = Covariation.collapse_sequences(id2seq_dict, refSeq, max_identity=0.98, min_match_identity=0.2, max_indel_ratio=0.5)
        clean_id2seq = Covariation.remove_bpbreak_sequences(uniq_seq, refStr, maxBpBreakRatio=0.2, maxBpDeletRatio=0.5, maxBpBreakCount=99, maxBpDeletion=99)
        ## 去除其中有非ATCG碱基的序列
        for key in list(clean_id2seq.keys()):
            v = clean_id2seq[key].replace('~','-').replace(':','-').replace('.','-').upper().replace('U','T')
            if len(v) != v.count('A')+v.count('T')+v.count('C')+v.count('G')+v.count('-'):
                del clean_id2seq[key]
        pure_fasta = { k:v.replace('~','').replace(':','').replace('.','').replace('-','').upper().replace('U','T') for k,v in clean_id2seq.items() }
        collapsed_pure_fasta = Alignment.cd_hit_est(pure_fasta, identity=0.98, verbose=True)
        clean_id2seq = { k:clean_id2seq[k] for k in collapsed_pure_fasta }
        if len(clean_id2seq)==0:
            #### No homologous sequence to call covariation in the sequence db
            return False
        ### 去除gap区域
        returnvalues = Covariation.remove_gap_columns(clean_id2seq, refSeq=refSeq, refStr=refStr, minGapRatio=1.0)
        returnvalues['id2seq'] = { k:v for k,v in returnvalues['id2seq'].items() if v.count('A') } # 去除其中有N的序列
        GS_DE = Covariation.read_sto_DE(cmsearch_output_sto)
        small_seqdb = { k:v.replace('-','') for k,v in returnvalues['id2seq'].items() }
        small_seqdb['input'] = input_seq
        General.write_fasta(small_seqdb, join(workdir, 'cmalign_input.fa'), GS_DE)
        h = Covariation.cmalign(input_cm, join(workdir, 'cmalign_input.fa'), cmalign_output_sto, cpu=0, glocal=False, 
            outformat='Stockholm', mxsize=1028.0, verbose=True, showCMD=True, use_LSF=True, LSF_parameters={'cpu': 5})
        h.wait()
        #os.system(f"cmalign -o {cmalign_output_sto} --outformat Stockholm {input_cm} {join(workdir, 'cmalign_input.fa')}")
        #os.remove(join(workdir, 'cmalign_input.fa'))
    return cmalign_output_sto

def align_long_secondary_str(virus_seq, virus_dot, seqDBfn, workdir, max_process=10, blastFilter=None):
    from multiprocessing import Process, Pool
    
    ### 1. 准备paramter list
    virus_stemloop = find_stemloop(virus_dot, max_dist=300)
    parameter_list = []
    for ls,re in virus_stemloop:
        if 10<re-ls+1<30:
            inputseq = virus_seq[ls-1:re]
            inputdot = virus_dot[ls-1:re]
            title = f'{ls}-{re}'
            fullpath = join(workdir, title)
            if os.path.exists(join(fullpath, "cmalign_output_2.sto")):
                continue
            if not os.path.exists(fullpath):
                os.mkdir(fullpath)
            Covariation.dot2sto({'input':[inputseq,inputdot]}, title, join(fullpath, 'input.sto'))
            params = (join(fullpath, 'input.sto'), inputseq, seqDBfn, fullpath, 3, blastFilter)
            parameter_list.append(params)
    
    #print([d[3] for d in parameter_list])
    pool = Pool(processes=max_process)
    pool.starmap(iterate_alignstr, parameter_list)

def find_stemloop(dot, max_dist=300):
    stemloop = []
    stems = Structure.find_stem(dot, max_stem_gap=3, min_stem_len=5)
    stems.sort()
    for ls,le,rs,re in stems:
        if rs-le<=max_dist and re-ls>=max_dist:
            stemloop.append( [ls,re] )
        elif rs-le>max_dist:
            continue
        else:
            if len([ d for d in stemloop if d[0]<ls<re<d[1] ]) == 0:
                stemloop.append([ls,re])
    return stemloop


virus_seq, virus_dot = General.load_dot('/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/SARS2.dot')['NC_045512.2']
seqDBfn = '/Share2/home/zhangqf7/tmp/covariation/GenomicFastaResults.uniq.fa'
workdir = '/Share2/home/zhangqf7/tmp/covariation/SASR2_validation'

blastFilter = { 
    'full_seq': virus_seq, 
    'blastdb':"/Share2/home/zhangqf7/tmp/covariation/blastdb/coronavirus", 
    'flanking': 20, 
    'blast_identity': 0.4 
}


align_long_secondary_str(virus_seq, virus_dot, seqDBfn, workdir, max_process=20, blastFilter=blastFilter)



