import pandas as pd
import os
import subprocess
df=pd.read_csv("../qcc/risultati/aggregate_sabre.csv")

for tm in [10,20,30,60]:
    for i in range(160):
        riga=df.iloc[i,:]
        sz=riga['size']
        if sz<100000:
            procs=[]
            if sz<500:
                nch=1
            elif sz<1000:
                nch=10
            elif sz<10000:
                nch=20
            elif sz<20000:
                nch=50
            elif sz<50000:
                nch=100
            else:
                nch=150
            #nch=1
            for s in range(1,11):
                cmd=f"./compile_qasm ../revlib/architectures/ibmq_tokyo.arch ../revlib/examples/{riga.instance_name}.qasm -timeout {tm} -bandit 1 -min_pr 0 -max_pr 0 -objf depth -divide {nch} -res esp_depth_bis.csv -seed {s}"
                elenco=cmd.split(' ')
                proc=subprocess.Popen(elenco)
                procs.append(proc)
            for p in procs:
                p.wait()
