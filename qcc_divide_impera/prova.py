import subprocess
lista_nch=[1,    2,    4,    5,   10,   20,   50,   75,  100,  200,  250, 500,  750, 1000]
circuit="ising_model_10"
tm=10
res_file="prova_sw.csv"
for nch in lista_nch:
    for s in range(1,15+1):
        cmd=f"./compile_qasm ../revlib/architectures/ibmq_tokyo.arch ../revlib/examples/{circuit}.qasm -timeout {tm} -bandit 1 -divide {nch} -res {res_file} -objf swaps -seed {s}"
        elenco=cmd.split(' ')
        proc=subprocess.Popen(elenco)
        proc.wait()

