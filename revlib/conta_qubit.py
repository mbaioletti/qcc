from qiskit import QuantumCircuit
import os

def conta_qubit(qc):
	used=set([])
	for d in qc.data:
		nq=d.operation.num_qubits
		for i in range(nq):
			j=d.qubits[i]._index
			used.add(j)
	return used
	

def scrivi_dati(name, file_output):
    qc=QuantumCircuit.from_qasm_file("examples/"+name)
    s=conta_qubit(qc)
    print(name, qc.size(), qc.depth(), len(s), file=file_output)
    
_,_,fs=next(os.walk("examples/"))

file_output=open("dati2.csv","w")
for f in fs:
    scrivi_dati(f, file_output)
    
file_output.close()
