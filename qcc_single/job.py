for s in `seq 1 10`
do 
./compile_qasm ../revlib/architectures/ibmq_tokyo.arch ../revlib/examples/qft_16.qasm  \
-timeout 60 -gamma 3 -beta 5 -delta 1  -min_pr 0.5 -max_pr 0.5 -choose rw_add \
-check 1 -seed $s -res prova.csv & 
done
