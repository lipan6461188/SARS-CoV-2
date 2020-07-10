
import General, Structure, Visual, SARS2

wuhan_id = 'NC_045512.2'
sequence = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')[wuhan_id]
shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-w50.shape')[wuhan_id]

def read_log10_Pb(pbfn):
    pairingProb = []
    for lc,line in enumerate(open(pbfn)):
        if lc>1:
            nc1,nc2,log10Prob = line.strip().split()
            pairingProb.append( (int(nc1),int(nc2),10**(-float(log10Prob))) )
    return pairingProb

def ctPlus(ct, plus):
    new_ct = []
    for bp in ct:
        if len(bp)==3:
            new_ct.append([bp[0]+plus, bp[1]+plus, bp[2]])
        else:
            new_ct.append([bp[0]+plus, bp[1]+plus])
    return new_ct

def refold_single_region(sequence, shape, dot, min_single=50):
    single_dot = "."*min_single
    pos = dot.find(single_dot)
    refold_bpprob = []
    while pos!=-1:
        end = pos
        while dot[end]=='.':
            end += 1
        print(f"refold {pos}-{end}")
        subseq = sequence[pos:end]
        subshape = shape[pos:end]
        subdot, subbpprob = SARS2.maxexpect_predict(subseq, subshape, si=-0.4, sm=1.5, md=400, remove_lowbp=True)
        #subdot = Structure.predict_structure(subseq, subshape, si=-0.4, sm=1.5)
        dot = list(dot)
        dot = dot[:pos]+list(subdot)+dot[end:]
        dot = "".join(dot)
        refold_bpprob += ctPlus(subbpprob, pos)
        pos = dot.find(single_dot)
    return dot, refold_bpprob

########################
## 1. Read pairing probability and structure
########################

root = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/full_length/'
filePrefix = ['1-5000', '4001-9000', '8001-13000', '12001-17000', '16001-21000', '20001-25000', '24001-29000', '24904-29903']
raw_ct_dict = {}
scaled_ct_dict = {}
log10Pb_dict = {}
for fp in filePrefix:
    start, end = fp.split('-')
    start = int(start)
    end = int(end)
    seq, ct, length = General.load_ct(join(root, f"{fp}.ct"))
    assert seq==sequence[start-1:end]
    assert length==5000
    
    tmp_dot = Structure.ct2dot(ct, 5000)
    tmp_prob = read_log10_Pb(join(root, f"{fp}.txt"))
    tmp_dot = SARS2.remove_lowporb_bp(tmp_prob, tmp_dot, minprob=0.6, minDist=150)
    bpprob = SARS2.collect_bpprob(tmp_prob, tmp_dot)
    
    #### Save ct
    raw_ct_dict[fp] = copy.deepcopy(ct)
    for i in range(len(bpprob)):
        bpprob[i] = (bpprob[i][0]+start-1, bpprob[i][1]+start-1, bpprob[i][2])
    scaled_ct_dict[fp] = bpprob

########################
## 2. Collect structure segment
########################

leftmost_base = 1
rightmost_base = None
full_ct_bpprob = []

last_dot = Structure.ct2dot(raw_ct_dict[filePrefix[0]], 5000)
last_sl_list = SARS2.find_stemloop(last_dot)

for ii in range(1, len(filePrefix)):
    
    cur_key = filePrefix[ii]
    cur_dot = Structure.ct2dot(raw_ct_dict[cur_key], 5000)
    cur_sl_list = SARS2.find_stemloop(cur_dot)
    plus = int(cur_key.split('-')[0])-1
    cur_sl_list = ctPlus(cur_sl_list, plus)
    
    i = len(last_sl_list)-1
    j = 0
    last_end = int(filePrefix[ii-1].split('-')[1])
    while cur_sl_list[j][0]<=last_end:
        j += 1
    
    j -= 1
    while last_sl_list[i][0]<=cur_sl_list[j][1] and cur_sl_list[j][0]<=last_sl_list[i][1]:
        i -= 1
    
    while j>=0 and last_sl_list[i][1]<cur_sl_list[j][0]:
        j -= 1
    
    print(i, j+1, last_sl_list[i], cur_sl_list[j+1])
    
    rightmost_base = last_sl_list[i][1]
    full_ct_bpprob += [ bp for bp in scaled_ct_dict[filePrefix[ii-1]] if leftmost_base<=bp[0]<bp[1]<=rightmost_base ]
    leftmost_base = cur_sl_list[j+1][0]
    
    last_sl_list = cur_sl_list

full_ct_bpprob += [ bp for bp in scaled_ct_dict[filePrefix[-1]] if leftmost_base<=bp[0] ]

########################
## 3. Refold and organize pairing probability
########################

full_dot = Structure.ct2dot([d[:2] for d in full_ct_bpprob], len(sequence))
refold_dot, refold_bpprob = refold_single_region(sequence, shape, full_dot, min_single=50)

### Check
based_bases = [ d[0] for d in full_ct_bpprob ] + [ d[1] for d in full_ct_bpprob ]
for b1,b2,v in refold_bpprob:
    assert b1 not in based_bases
    assert b2 not in based_bases

full_ct_bpprob += refold_bpprob
full_ct_bpprob.sort(key=lambda x: x[0])

### Check each base pair have bp prob
tmp_ct = Structure.dot2ct(refold_dot)
i = 0
assert len(tmp_ct)==len(full_ct_bpprob)
while i<len(tmp_ct):
    assert list(tmp_ct[i][:2])==list(full_ct_bpprob[i][:2]), f"Error {i} {tmp_ct[i]} {full_ct_bpprob[i]}"
    i += 1

########################
## 4. Save dot and pairing probability
########################

structure_fn = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/SARS2.dot'
bpprob_fn = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/SARS2_bpprob.txt'

General.write_dot({'NC_045512.2':[sequence,refold_dot]}, structure_fn)
print("\n".join([f"{b1}\t{b2}\t{v:.3f}" for b1,b2,v in full_ct_bpprob]), file=open(bpprob_fn, 'w'))

