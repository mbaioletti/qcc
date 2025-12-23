from ctypes import cdll, c_int, c_wchar_p, Structure, POINTER, byref, pointer
from qiskit import QuantumCircuit
from qiskit.circuit import library
lib=cdll.LoadLibrary("./libcirc.so")
import numpy as np

class DirshGate(Structure):
     _fields_ = [("name", c_wchar_p),
                 ("num_qubits", c_int),
                 ("qb1", c_int),
                 ("qb2", c_int),
                 ("param", c_wchar_p)
                ]
                
lib.dirsh.restype = POINTER(DirshGate)

def convert_to_dirsh(qc):
    n=len(qc.data)
    gates_creator=DirshGate*n
    c=gates_creator()
    for i in range(n):
        c[i].name=qc.data[i].operation.name
        c[i].num_qubits=qc.data[i].operation.num_qubits
        c[i].qb1=qc.data[i].qubits[0]._index
        c[i].qb2=qc.data[i].qubits[1]._index if qc.data[i].operation.num_qubits==2 else -1
        params=qc.data[i].operation.params
        if params != []:
            c[i].param=str(params[0])
        else:
            c[i].param=""
    return c

def find_num_qubits(c1):
    n=0
    i=0
    while c1[i].name != "end":
        n=max(n,c1[i].qb1)
        if c1[i].num_qubits==2:
            n=max(n,c1[i].qb2)
        i +=1
    return n+1

gate_list={
    "x": library.XGate,
    "y": library.YGate,
    "z": library.ZGate,
    "h": library.HGate,
    "t": library.TGate,    
    "s": library.SGate,        
    "cx": library.CXGate,        
    "swap": library.SwapGate,
    "rx": library.RXGate,            
    "ry": library.RYGate,            
    "rz": library.RZGate                    
}
    
def convert_to_qiskit(c1):
    nqb=find_num_qubits(c1)
    qc2=QuantumCircuit(nqb)
    i=0
    while c1[i].name != "end":
        g=c1[i].name.lower()
        if g in gate_list:
            gg=gate_list[g]
            l=[c1[i].qb1]
            if c1[i].num_qubits==2 or g=="swap":
                l.append(c1[i].qb2)
            if c1[i].param != "":
                qc2.append(gg(float(c1[i].params)), l)
            else:
                qc2.append(gg(), l)
        else:
            print(f"Gate {g} not found")
            exit()
        i += 1
    return qc2
    
def run_dirsh(qc, fname_arch):
    c=convert_to_dirsh(qc)
    c1=lib.dirsh(c, len(qc.data), fname_arch)
    qc2=convert_to_qiskit(c1)
    return qc2
    

qc=QuantumCircuit(2,2)
qc.x(0)
qc.h(1)
qc.cx(0,1)

qc2=run_dirsh(qc, "../revlib/architectures/ibmq_tokyo.arch")
