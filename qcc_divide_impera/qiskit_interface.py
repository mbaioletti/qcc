from ctypes import cdll, c_int, c_wchar_p, Structure, POINTER, byref, pointer
from qiskit import QuantumCircuit
lib=cdll.LoadLibrary("./libcirc.so")
import numpy as np

class DirshGate(Structure):
     _fields_ = [("name", c_wchar_p),
                 ("num_qubits", c_int),
                 ("qb1", c_int),
                 ("qb2", c_int)
                ]
                
lib.dirsh.restype = POINTER(DirshGate)

def run_dirsh(qc, fname_arch):
    n=len(qc.data)
    gates_creator=DirshGate*n
    c=gates_creator()
    for i in range(n):
        c[i].name=qc.data[i].operation.name
        c[i].num_qubits=qc.data[i].operation.num_qubits
        c[i].qb1=qc.data[i].qubits[0]._index
        c[i].qb2=qc.data[i].qubits[1]._index if qc.data[i].operation.num_qubits==2 else -1
    c1=lib.dirsh(c, n, fname_arch)
    return c1
    

qc=QuantumCircuit(2,2)
qc.x(0)
qc.h(1)
qc.cx(0,1)

qc2=run_dirsh(qc, "../revlib/architectures/ibmq_tokyo.arch")
