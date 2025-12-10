import pandas as pd
import os
import subprocess
df=pd.read_csv("../revlib/tutti.csv", sep=",")

files=["miller_11.qasm", "4gt13_91.qasm", "0410184_169.qasm", "ex3_229.qasm", "sf_276.qasm", "rd53_131.qasm","rd53_251.qasm", "sqrt8_260.qasm", "misex1_241.qasm"]

for tm in [10]:
    for f in files:
        riga=df[df.instance_name==f].iloc[0,:]
        sz=int(riga['size'])
        de=int(riga.depth)
        qb=int(riga['qubits'])
        if sz<100000:
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
            print(f"{f}, {sz}, {qb}")
            for si in range(0,50+1):
                procs=[]
                for sa in range(1,15+1):
                    cmd=f"./compile_qasm ../revlib/architectures/ibmq_tokyo.arch ../revlib/examples/{f} -timeout {tm} -bandit 1 -divide {nch} -res stato_iniz.csv -seed {sa} -seed_is {si}"
                    elenco=cmd.split(' ')
                    proc=subprocess.Popen(elenco)
                    procs.append(proc)
                for p in procs:
                    p.wait()

