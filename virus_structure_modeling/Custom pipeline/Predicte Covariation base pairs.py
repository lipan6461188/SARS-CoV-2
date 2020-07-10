
import General, Structure, Visual, SARS2

wuhan_id = 'NC_045512.2'
sequence = General.load_fasta('/Share2/home/zhangqf7/lipan/SARS2/sequence/SARS2.fa')[wuhan_id]
shape = General.load_shape('/Share2/home/zhangqf7/lipan/SARS2/icSHAPE/2020-06-01-process/virus-w50.shape')[wuhan_id]
dot = General.load_dot('/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/SARS2.dot')[wuhan_id][1]

seqdbFn = "/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/homoseq/Rfam_sequence.fa"
workdir_root = "/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/covariation"
covary_bases = SARS2.call_covariation_split(sequence, dot, seqdbFn, workdir_root)




