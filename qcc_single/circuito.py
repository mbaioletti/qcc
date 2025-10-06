from ctypes import cdll, c_int, c_wchar_p, Structure
from qiskit import QuantumCircuit
lib=cdll.LoadLibrary("./libcirc.so")
import numpy as np

class Coppia(Structure):
     _fields_ = [("name", c_wchar_p),
                 ("num_qubits", c_int),
                 ("qb1", c_int),
                 ("qb2", c_int)
                ]

def invia_circuito(qc):
    n=len(qc.data)
    coppie=Coppia*n
    c=coppie()
    for i in range(n):
        c[i].name=qc.data[i].operation.name
        c[i].num_qubits=qc.data[i].operation.num_qubits
        c[i].qb1=qc.data[i].qubits[0]._index
        c[i].qb2=qc.data[i].qubits[1]._index if qc.data[i].operation.num_qubits==2 else -1
    lib.leggi_circuito(c, n)

qc=QuantumCircuit.from_qasm_file("soluzione_sabre.qasm")
invia_circuito(qc)

