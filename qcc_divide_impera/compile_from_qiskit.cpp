#include "compiler.h"
#include <cwchar>

struct Gate_from_qiskit {
	const wchar_t *name;
	int num_qubits;
	int qb1, qb2;
};

string from_wchar(const wchar_t *c) {
	string s="";
	for(int j=0; c[j]!=0; j++)
		s += char(c[j]);
	return s;
}

Gate_from_qiskit* call_dirsh(Gate_from_qiskit c[], int n, string arch_file) {
    Problem p;
    p.num_gates=0;
    int last_qubit=-1;
	for(int i=0; i<n; i++) {
		string type=from_wchar(c[i].name);
        int qb1=c[i].qb1, qb2=c[i].qb2;
        string param="";
        p.gates.push_back(Gate(p.num_gates,type,qb1,qb2,param));
        last_qubit=max(last_qubit, qb1);
        last_qubit=max(last_qubit, qb2);        
        p.num_gates++;
	}
    p.n_log_qubits = last_qubit+1;
    p.compute_precedences();
    cout << "Read circuit  with " << p.gates.size() << " gates with depth " << p.depth << endl;
    p.load_architecture(arch_file);
    compile_data.problem = &p;
    compile_data.timeout=1;
    Solution *s=optimize();
    if(s==nullptr) {
        cout << "No solution found" << endl;
    }
    else {
        cout << "Found a solution with " << s->activities.size() << " gates, depth " << s->makespan << " and " << s->num_swaps << " swaps " << endl;
    }
    int ng=s->activities.size();
    Gate_from_qiskit *res=new Gate_from_qiskit[ng+1];
    for(int i=0; i<ng; i++) {
        wstring w(s->activities[i]->gate->type.begin(), s->activities[i]->gate->type.end());
        res[i].name=w.c_str();
    }
    
    return res;
}

extern "C" {
	Gate_from_qiskit* dirsh(Gate_from_qiskit c[], int n, wchar_t *fa) {
		return call_dirsh(c, n, from_wchar(fa));
	}
}
