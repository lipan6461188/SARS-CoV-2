
splice_sites = """
64  28254   7012631 28273
65  21551   1117576 21562
66  27384   688112  27393
65  25380   618367  25392
65  27883   229749  27893
64  26467   215361  26522
69  26236   149744  26244
69  27040   38086   27201
68  28262   37740   28273
63  21549   17894   21562
63  28253   13079   28273
76  28266   9868    28273
60  27378   9320    27393
71  27761   8568    27824
61  28250   6457    28273
68  15776   6103    15811
70  27760   5982    69
69  27673   3866    27755
63  27881   3702    27893
61  28251   2946    28273
65  27383   2814    27393
66  27484   1901    27755
68  29153   1216    29221
76  21562   992 21562
74  21057   990 21070
76  27394   950 27755
83  21562   863 21562
74  21552   843 21562
68  29111   818 29221
67  24776   728 24809
70  22944   703 23003
70  22276   699 22358
68  28282   686 28283
82  26276   669 26522
75  28964   622 28972
67  26865   606 26895
70  22501   580 22613
71  27675   573 27755
62  28254   554 28273
69  24890   553 24890
70  29002   514 29078
68  27676   492 27755
71  26505   489 26522
68  21549   427 21562
75  26485   409 26522
63  25378   408 25392
70  21305   396 21316
69  18940   360 18940
68  27043   360 27201
78  25393   355 25404
69  29165   345 29221
68  27762   338 27824
72  21546   290 21562
64  26497   268 26522
63  5784    264 5858
69  23027   263 23045
68  25424   258 25456
77  26271   258 26522
69  19570   256 19630
75  27798   255 27824
57  21539   253 21562
70  20315   250 20329
63  26466   230 26522
68  26235   230 26244
64  27382   225 27393
67  2291    222 2455
70  14371   218 14380
70  2793    217 69
77  26480   212 26522
73  23027   210 23045
66  27298   204 27358
68  26289   188 26522
59  28249   186 28273
69  19111   185 19147
67  22404   185 22406
69  28164   184 28206
83  26250   183 26522
61  26464   179 26522
68  5788    178 5858
68  18050   175 18152
74  28265   175 28273
68  28994   158 29078
68  27866   153 27867
69  28206   151 28206
67  29378   149 29402
72  28570   146 28573
75  28267   144 28273
69  22558   144 22613
67  22502   143 22613
71  24892   137 24941
78  27754   136 27755
69  27758   135 27824
73  27777   134 27824
67  26925   134 27201
68  23364   133 23402
68  25490   133 25523
68  27473   133 27755
66  23849   131 24050
74  28891   128 28900
69  29590   126 29617
""".strip().split('\n')
splice_sites = [ d.split() for d in splice_sites ]
splice_sites.sort(key=lambda x: int(x[2]), reverse=True)
splice_rank = { data[0]+"_"+data[1]:rank+1 for rank,data in enumerate(splice_sites) }

splice_seq = General.load_fasta('/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/14.map_virus_splice-all/index/splice_sequence.fa')
splice_shape_fn = '/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/14.map_virus_splice-all/only_cover_splice/final.shape'
splice_shape = General.load_shape(splice_shape_fn)

Dot = {}
for key in splice_shape:
    if splice_shape[key][:300].count('NULL')<30:
        subseq = splice_seq[key][:300]
        subshape = splice_shape[key][:300]
        Dot[key] = Structure.predict_structure(subseq, subshape)

data = []
for key in Dot:
    subseq = splice_seq[key][:300]
    subshape = splice_shape[key][:300]
    ss_num = Dot[key][:100].count('.')
    data.append([ss_num, key, subseq, Dot[key], subshape])

data.sort(key=lambda x:x[0], reverse=True)

OUT = open(join(HOME,"figs/splice_structure.ps1"), 'w')
rank = 1
for ss_num, key, subseq, dot, subshape in data:
    AUC = General.calc_AUC_v2(dot, subshape)
    title = f"{key} str_rank={rank} splice_rank={splice_rank[key]} AUC={AUC:.3}"
    cmd = Visual.Plot_RNAStructure_Shape(subseq, dot, subshape, title=title, mode='label')
    print(cmd, file=OUT)
    rank += 1

OUT.close()





