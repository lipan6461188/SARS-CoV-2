
import General, Structure, Visual, SARS2

wuhan_id = 'NC_045512.2'

shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-w50.shape')[wuhan_id]
sequence = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')[wuhan_id]

Len = len(sequence)
out = '/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/full_length'
start = 0
finish = False
winsize = 5000
overlap = 1000
while not finish:
    if start+winsize>Len:
        start = Len - winsize
        finish = True
    shapefn = f"{start+1}-{start+winsize}.shape"
    seqfn = f"{start+1}-{start+winsize}.seq"
    Structure.__build_SHAPE_constraint(shape[start:start+winsize], join(out, shapefn))
    Structure.__build_single_seq_fasta(sequence[start:start+winsize], join(out, seqfn))
    start += winsize - overlap

