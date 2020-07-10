
##### 1. Partition

filePrefix=(1-5000 4001-9000 8001-13000 12001-17000 16001-21000 20001-25000 24001-29000 24904-29903)
for fp in ${filePrefix[@]};
do
    bsub -q Z-ZQF -e error -n 20 \
        "partition-smp ${fp}.seq \
            ${fp}.pfs \
            -si -0.4 -sm 1.5 \
            --SHAPE ${fp}.shape \
            --maxdistance 400"  
done

##### 2. ProbabilityPlot

for fp in ${filePrefix[@]};
do
    bsub -q Z-ZQF -e error -n 20 \
        "ProbabilityPlot ${fp}.pfs ${fp}.txt --text"  
done

##### 3. MaxExpect

for fp in ${filePrefix[@]};
do
    bsub -q Z-ZQF -e error -n 20 \
        "MaxExpect-smp ${fp}.pfs ${fp}.ct --structures 10"  
done






