
import General, Structure, Visual, SARS2, Covariation
importCommon()

####################
### Read Sequence, Shape, Dot-bracket structure
####################

wuhan_id = 'NC_045512.2'
sequence = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')[wuhan_id]
shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-w50.shape')[wuhan_id]
dot = General.load_dot('/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/SARS2.dot')[wuhan_id][1]

####################
### Read base pairing probability
####################

bpfn = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/SARS2_bpprob.txt'
bpprob = []
for line in open(bpfn):
    data = line.strip().split()
    b1, b2, v = int(data[0]), int(data[1]), float(data[2])
    bpprob.append( (b1, b2, v) )

####################
### Read Covariation
####################

def batch_read_cov_RScape(virus_sequence):
    infolder = "/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/covariation-without-nohmm"
    folders = os.listdir(infolder)
    full_covary_bps = []
    for folder in folders:
        combined_folder = join(infolder, folder)
        if os.path.isdir(combined_folder):
            start, end = folder.split('-')
            start, end = int(start), int(end)
            query_seq = virus_sequence[start-1:end]
            
            rscape_file = os.path.join(combined_folder, 'R-scape', f'{folder}.cov')
            covary_bps = []
            if os.path.exists(rscape_file):
                rscape_list = Covariation.read_RScape_result(rscape_file)
                id2seq, refStr, refAnnot = General.load_stockholm( join(combined_folder, 'output.sto') )[0]
                input_id = [ key for key in id2seq if key.startswith('input') ][0]
                posDict = Covariation.get_alignedPos2cleanPos_dict(id2seq[input_id])
                for bp in rscape_list:
                    covary_bps.append( [posDict[bp[0]], posDict[bp[1]]] )
                    l, r = covary_bps[-1]
                    bp_base = query_seq[l-1]+query_seq[r-1]
                    if bp_base.replace('U','T') not in ('AT','TA','GC','CG','GT','TG'):
                        print(f"Warning: {bp_base} is not a RNA base pair")
            else:
                print(f"Warning: {rscape_file} doesn't exists")
            for i in range(len(covary_bps)):
                covary_bps[i][0] += start-1
                covary_bps[i][1] += start-1
            full_covary_bps += covary_bps
    full_covary_bps.sort()
    for bp in full_covary_bps:
        bp_base = sequence[bp[0]-1]+sequence[bp[1]-1]
        assert bp_base in ('AT','TA','GC','CG','GT','TG')
    print(f"Size of covarying base pairs: {len(full_covary_bps)}")
    return full_covary_bps

#full_covary_bps = batch_read_cov_RScape(sequence)

def batch_read_cov_RNAalignfold(virus_sequence):
    infolder = "/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/covariation-without-nohmm"
    folders = os.listdir(infolder)
    full_covary_bps = []
    bar = tqdm(total=len(folders), leave=True)
    for folder in folders:
        bar.update(1)
        combined_folder = join(infolder, folder)
        if os.path.isdir(combined_folder):
            start, end = folder.split('-')
            start, end = int(start), int(end)
            query_seq = virus_sequence[start-1:end]
            
            clean_sto_file = os.path.join(combined_folder, 'R-scape', f'{folder}.cons.sto')
            covary_bps = []
            if os.path.exists(clean_sto_file):
                covary_bps = Covariation.calc_covBP_from_sto(clean_sto_file, query_seq, min_score=0.4)
                for i in range(len(covary_bps)):
                    covary_bps[i][0] += start-1
                    covary_bps[i][1] += start-1
                full_covary_bps += covary_bps
    bar.close()
    full_covary_bps.sort()
    for bp in full_covary_bps:
        bp_base = sequence[bp[0]-1]+sequence[bp[1]-1]
        assert bp_base in ('AT','TA','GC','CG','GT','TG'), f"{bp}: {bp_base}"
    print(f"Size of covarying base pairs: {len(full_covary_bps)}")
    return full_covary_bps

full_covary_bps = batch_read_cov_RNAalignfold(sequence)


####################
### Generate CMD command
####################

stemloops = SARS2.find_stemloop(dot)

def is_large_sl(sl):
    if sl[1]-sl[0]>200:
        return True
    return False

def manual_period(raw_cmd, seqlen, cur_pos, priod):
    cmd = raw_cmd + " -periodNum 0 "
    j = 0
    annotation_cmd = ""
    for i in range(cur_pos, cur_pos+seqlen):
        j += 1
        if i%priod == 0:
            annotation_cmd += f"{i}:type=B,anchor={j},size=8,color=#ff5722;"
    if annotation_cmd:
        cmd += f"-annotations \"{annotation_cmd}\""
    return cmd

start = 1
i = 0
OUT = open(join(HOME, 'figs/plot_str-fill.ps1'), 'w')
while i<len(stemloops):
    while i<len(stemloops)-1 and not is_large_sl(stemloops[i+1]) and stemloops[i+1][1]-start<1000:
        i += 1
    end = stemloops[i][1]
    
    ### Prepare Covary base pairs
    cur_covarying_bps = [ bp for bp in full_covary_bps if start<=bp[0]<bp[1]<=end ]
    highlight_region = []
    for bp in cur_covarying_bps:
        if bp[2]<0.5:
            color = "#08b0f2"
        elif bp[2]<0.7:
            color = "#12c634"
        else:
            color = "#ed0c0c"
        
        highlight_region.append( [bp[0]-start+1, bp[0]-start+1, color] )
        highlight_region.append( [bp[1]-start+1, bp[1]-start+1, color] )
    
    subbpprob = [ (d[0]-start+1,d[1]-start+1,d[2]) for d in bpprob if start<=d[0]<d[1]<=end ]
    title = f"{start}-{end}"
    cmd = Visual.Plot_RNAStructure_Shape(sequence[start-1:end], 
        dot[start-1:end],
        shape[start-1:end],
        mode='fill',
        highlight_region=highlight_region,
        bpprob=subbpprob,
        bpprob_mode='both',
        title=title,
        bpwarning=False,
        period=100, first_base_pos=start, peroid_color='#ff5722')
    #cmd = manual_period(cmd, end-start+1, start, 100)
    print(cmd, file=OUT)
    start = end+1
    i += 1

OUT.close()

####################
### Frame-Shift Element
####################

subseq = sequence[13476-1:13503]
subshape = shape[13476-1:13503]
subdot = Structure.predict_structure(subseq)
cmd = Visual.Plot_RNAStructure_Shape(subseq, subdot, subshape, mode='fill', first_base_pos=13476)
print(cmd)

