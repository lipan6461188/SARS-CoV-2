
cd /Share2/home/zhangqf7/lipan/SARS2/SM_folding_pipeline

### 1. 生成预测命令

bash runManyParts.sh SARS2 > commands.sh

### 2. 使用partition预测配对概率

bash commands.sh

### 3. 合并ShannoEntropy

cp recombineShannon.sh PFS_Final/
cd PFS_Final/
recombineShannon.sh 29903

### 4. 收集大于99%配对概率的碱基对做成ELEMENTS.ct

python masterModel.py SARS2.ct PFS_Final/ ELEMENTS.ct

### 5. 使用前面的结构作为碱基对限制使用Fold预测新的二级结构

cp RNAFILE.shape ELEMENTS.shape
bash runFinalFold.sh ELEMENTS > commands2.sh
bash commands2.sh

### 6. 合并分别预测的文件

python masterModel.py ELEMENTS.ct FinalFold/ finalStructure.ct

