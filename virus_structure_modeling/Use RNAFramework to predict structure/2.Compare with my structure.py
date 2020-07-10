

shape = General.load_shape('/Share/home/zhangqf8/sunlei/data/SARS2/20200529-SARS2-huh7/SARS2-huh7/Processing/8.calcSHAPE/virus-w50.shape')['NC_045512.2']

seq,dot_rf = General.load_dot("/Share2/home/zhangqf7/lipan/SARS2/RNAFramework/rf_fold2/structures/SARS2.db")['SARS2']
#seq,ctList,length = General.load_ct("/Share2/home/zhangqf7/lipan/SARS2/SM_folding_pipeline/finalStructure.ct")
#dot_SM = Structure.ct2dot(ctList, length)

### 先看看5'UTR

cmd = Visual.Plot_RNAStructure_Shape(seq[:300], dot_rf[:300], shape[:300])
print(cmd)

### 再看看3'UTR

cmd = Visual.Plot_RNAStructure_Shape(seq[29615:], dot_rf[29615:], shape[29615:])
print(cmd)

### 和我预测的结构做比较

sars2_dot = General.load_dot("/Share2/home/zhangqf7/lipan/SARS2/predict_structure/SARS2/SARS2.dot")['NC_045512.2'][1]

Structure.evaluate_dot(sars2_dot, dot_rf)

# 有多少个公共的碱基对

ct1 = Structure.dot2ct(dot_rf)
ct2 = Structure.dot2ct(sars2_dot)

print(len(ct1))
len(set(ct1) & set(ct2))
print(len(ct2))


