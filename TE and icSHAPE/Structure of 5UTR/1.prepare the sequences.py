

import General, Structure, Visual

wuhan_id = 'NC_045512.2'
#shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus.shape')[wuhan_id]
seq = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')[wuhan_id]

link_locations = """
64  28254
65  21551
66  27384
65  25380
65  27883
64  26467
69  26236
69  27040
68  28262
63  21549
""".strip().split('\n')
link_locations = [ [int(d.split()[0]), int(d.split()[1])] for d in link_locations ]

ref_sequences = {}
for p1,p2 in link_locations:
    ref_sequences[ f"{p1}_{p2}" ] = seq[:p1]+seq[p2:p2+300]

inFn = "/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/splice_fasta/splice_sequence.fa"
General.write_fasta(ref_sequences, inFn)






