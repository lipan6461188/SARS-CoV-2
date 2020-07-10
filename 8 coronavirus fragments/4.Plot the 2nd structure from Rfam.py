
utr_sequences = General.load_fasta("/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/sequences/utr_sequences.fa")

#########################
### Code to read cmscan result
########################

class CmscanHit:
    def __init__(self):
        self.b_NC = None
        self.b_str = None
        self.b_ref = None
        self.b_score = None
        self.b_query = None
        self.b_PP = None
        
        self.ref_range = None
        self.query_range = None
        self.ref_name = None
        self.query_name = None
    
    def check(self):
        assert len(self.b_NC)==len(self.b_str), "len(self.b_NC)!=len(self.b_str)"
        assert len(self.b_NC)==len(self.b_ref), "len(self.b_NC)!=len(self.b_ref)"
        assert len(self.b_NC)==len(self.b_score), "len(self.b_NC)!=len(self.b_score)"
        assert len(self.b_NC)==len(self.b_query), "len(self.b_NC)!=len(self.b_query)"
        assert len(self.b_NC)==len(self.b_PP), "len(self.b_NC)!=len(self.b_PP)"
    
    def __repr__(self):
        return f"{self.query_name}:{self.query_range[0]}-{self.query_range[1]} <==> {self.ref_name}:{self.ref_range[0]}-{self.ref_range[1]}"

def read_a_hit(IN):
    """
    Read cmscan -o file
    """
    line = IN.readline()
    while line and (not line.startswith(">>")):
        line = IN.readline()
    if not line:
        return False
    #assert line.startswith('>>'), "read_a_hit error: Hit should startswith >>"
    ref_name = line.strip()
    
    ### Here Only consider one hit
    for _ in range(4):
        IN.readline()
    NC_line = IN.readline().rstrip()
    str_line = IN.readline().rstrip()
    ref_line = IN.readline().rstrip()
    score_line = IN.readline().rstrip('\n')
    query_line = IN.readline().rstrip()
    PP_line = IN.readline().rstrip()
    #IN.readline()
    
    ## Count the blank
    bc = 0
    while PP_line[bc] == " ":
        bc += 1
    
    hit = CmscanHit()
    hit.b_NC = NC_line[bc:-3]
    hit.b_str = str_line.strip().split()[0]
    
    ref_name, ref_start, ref_seq, ref_end = ref_line.strip().split()
    hit.b_ref = ref_seq
    hit.ref_range = [ int(ref_start), int(ref_end) ]
    hit.ref_name = ref_name
    
    query_name, query_start, query_seq, query_end = query_line.strip().split()
    hit.b_query = query_seq
    hit.query_range = [ int(query_start), int(query_end) ]
    hit.query_name = query_name
    
    hit.b_score = score_line[bc:]
    hit.b_PP = PP_line.strip().split()[0]
    
    hit.check()
    return hit


IN = open("/Share2/home/zhangqf7/lipan/SARS2/Target-icSHAPE-MaP/sequences/cmscan_result.out")
hits = []
hit = read_a_hit(IN)
while hit:
    hits.append(hit)
    hit = read_a_hit(IN)

#########################
### Remove gap from raw dot-bracket
########################

def getstr(hit):
    new_str = ""
    new_seq = ""
    b_str = hit.b_str.replace(',','.').replace(':','.').replace('_', '.').replace('-', '.')
    for nc,base,d in zip(hit.b_NC,hit.b_query.upper(),b_str):
        if base != '-':
            if nc == 'v':
                new_str += '.'
            else:
                new_str += d
            new_seq += base
    new_str = Structure.ct2dot(Structure.dot2ct(new_str),len(new_str))
    return new_seq, new_str

results = []
for hit in hits:
    c1 = hit.query_name in ('MERS','bHKU9','hHKU1','bHKU5') and hit.ref_name in ('bCoV-5UTR') #('Corona_pk3', 'bCoV-5UTR', 'bCoV-3UTR')
    c2 = hit.query_name in ('SARS2-C','SARS2-T','SARS') and hit.ref_name in ('Sarbecovirus-5UTR') #('Corona_pk3', 's2m', 'Sarbecovirus-5UTR', 'Sarbecovirus-3UTR')
    c3 = hit.query_name in ('hNL63') and hit.ref_name in ('aCoV-5UTR') #('Corona_pk3', 'aCoV-5UTR', 'aCoV-3UTR')
    if c1 or c2 or c3:
        q_seq, q_str = getstr(hit)
        results.append([hit.query_name, hit.ref_name, hit.query_range[0], hit.query_range[1], q_seq, q_str])

#########################
### Read icSHAPE-MaP
########################

def read_UTR_icSHAPE():
    root = "/Share/home/zhangqf8/sunlei/data/SARS2/20200623-inVitro-SRAS-huh7-8UTR/8UTR-icSHAPE"
    icSHAPE = {}
    icSHAPE.update( General.load_shape(join(root, "MERS-HKU9/OUT/9.transSHAPE/virus.shape")) )
    icSHAPE.update( General.load_shape(join(root, "SARS2-C-NL63_HKU1/OUT/9.transSHAPE/virus.shape")) )
    icSHAPE.update( General.load_shape(join(root, "SARS2-T/OUT/9.transSHAPE/virus.shape")) )
    icSHAPE.update( General.load_shape(join(root, "SARS-HKU5/OUT/9.transSHAPE/virus.shape")) )
    return icSHAPE

shape_score = read_UTR_icSHAPE()



#########################
### Plot and annotation
########################

refToRF = {
    'Corona_pk3':       'RF00165',
    's2m':              'RF00164',
    'Corona_package':   'RF00182',
    'Corona_FSE':       'RF00507',
    'aCoV-5UTR':        'RF03116',
    'bCoV-5UTR':        'RF03117',
    'gCoV-5UTR':        'RF03118',
    'dCoV-5UTR':        'RF03119',
    'Sarbecovirus-5UTR':'RF03120',
    'aCoV-3UTR':        'RF03121',
    'bCoV-3UTR':        'RF03122',
    'gCoV-3UTR':        'RF03123',
    'dCoV-3UTR':        'RF03124',
    'Sarbecovirus-3UTR':'RF03125'
}


OUT = open(join(HOME, "figs/icSHAPE-5UTR.ps1"), 'w')
for name,ref_name,start,end,subseq,subdot in results:
    assert subseq.replace('U','T')==utr_sequences[name][start-1:end]
    if name not in shape_score: continue
    shape = shape_score[name][start-1:end]
    auc = round(General.calc_AUC_v2(subdot, shape),3)
    title = f"{name} {start}-{end} ; {ref_name}({refToRF[ref_name]}); {auc}"
    cmd = Visual.Plot_RNAStructure_Shape(subseq, subdot, shape, title=title)
    fname = name+"-"+refToRF[ref_name]+".png"
    print(cmd, file=OUT)

OUT.close()



