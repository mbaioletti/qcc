import pandas as pd
import subprocess
import numpy as np

df=pd.read_csv("../revlib/dati2.csv")

def esegui(fname, init_state):
    procs=[]
    for s in range(1,11):
        cmd=f"./compile_qasm ../revlib/architectures/ibmq_tokyo.arch ../revlib/examples/{fname} -timeout 10 -min_pr 0 -max_pr 0 -objf depth -res esp_angelo_1.csv -seed {s} -init_state {init_state}"
        elenco=cmd.split(' ')
        proc=subprocess.Popen(elenco)
        procs.append(proc)
    for p in procs:
        p.wait()
    print("task completato")


fname="rd53_130.qasm"

nqubits=7

for orders in range(30):
    init_state=list(range(nqubits))

    np.random.shuffle(init_state)

    init_state=";".join(str(x) for x in init_state)

    esegui(fname, init_state)
