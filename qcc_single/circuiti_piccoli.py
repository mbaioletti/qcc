import pandas as pd
import os
import subprocess
df=pd.read_csv("risultati/confronto_depth_sabre_rw_add.csv")

for i in range(160):
    riga=df.iloc[i,:]
    if riga['size']<6000:
        procs=[]
        for s in range(1,11):
            cmd=f"./compile_qasm ../revlib/architectures/ibmq_tokyo.arch ../revlib/examples/{riga.qasm}.qasm -timeout 60 -beta 0 -gamma 3 -delta 1 -res circuiti_piccoli_swap.csv -seed {s}"
            elenco=cmd.split(' ')
            proc=subprocess.Popen(elenco)
            procs.append(proc)
        for p in procs:
            p.wait()

