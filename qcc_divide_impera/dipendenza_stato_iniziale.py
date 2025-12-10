import pandas as pd
import os
import subprocess
df=pd.read_csv("../revlib/tutti.csv", sep=",")

files=["rd53_131","rd53_251"]

for tm in [10]:
    for f in files:
        riga=df[df.instance_name==(f+".qasm")]
        sz=int(riga['size'])
        de=int(riga.depth)
        qb=int(riga['qubits'])
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
            print(f"{f}, {sz}, {qb}")
"""            for s in range(1,11):
                cmd=f"./compile_qasm ../revlib/architectures/ibmq_tokyo.arch ../revlib/examples/{riga.instance_name}.qasm -timeout {tm} -bandit 1 -divide {nch} -res esp_swap_one.csv -seed {s}"
                elenco=cmd.split(' ')
                proc=subprocess.Popen(elenco)
                procs.append(proc)
            for p in procs:
                p.wait()
"""
