
cd /Share2/home/zhangqf7/lipan/SARS2/RNAFramework
#  这里去掉了-nlp，因为winfold时c容易报错
rf-fold \
    -o rf_fold \
    -ow \
    --folding-method 1 \
    --processors 20 \
    --slope 1.5 \
    --intercept -0.4 \
    -w -fw 3000 -fo 300 -wt 200 -pw 1500 -po 250 -dp -sh -md 600 \
    -KT \
    input.xml

