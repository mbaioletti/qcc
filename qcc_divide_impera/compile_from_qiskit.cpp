#include "compiler.h"

struct Gate_from_qiskit {
	wchar_t *name;
	int num_qubits;
	int qb1, qb2;
};

string from_wchar(wchar_t *c) {
	string s="";
	for(int j=0; c[j]!=0; j++)
		s += char(c[j]);
	return s;
}

void compile_orig(Gate_from_qiskit c[], int n) {
    Problem p;
    p.num_gates=0;
    int last_qubit=-1;
    cout << "Ci sono " << n << " gate sul circuito" << endl;
	for(int i=0; i<n; i++) {
		string type=from_wchar(c[i].name);
        cout << "Try to insert " << type << endl;
        int qb1=c[i].qb1, qb2=c[i].qb2;
        string param="";
        p.gates.push_back(Gate(p.num_gates,type,qb1,qb2,param));
        last_qubit=max(last_qubit, qb1);
        last_qubit=max(last_qubit, qb2);        
        p.num_gates++;
	}
    p.n_log_qubits = last_qubit+1;
    cout << "logical qubits " << p.n_log_qubits << endl;
    p.compute_precedences();
    cout << "Read circuit  with " << p.gates.size() << " gates with depth " << p.depth << endl;
}

extern "C" {
	void compile(Gate_from_qiskit c[], int n) {
		compile_orig(c, n);
	}
}
