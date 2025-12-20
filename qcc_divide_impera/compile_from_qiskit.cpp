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

wchar_t* to_wchar(const string &s) {
    int n=s.length();
    wchar_t* r=new wchar_t[n+1];
    for(int i=0; i<n; i++)
        r[i]=s[i];
    r[n]=0;
    return r;
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
    Gate_from_qiskit *res=new Gate_from_qiskit[1];
    res[0].name=to_wchar("prova");
    return res;
}

extern "C" {
	Gate_from_qiskit* dirsh(Gate_from_qiskit c[], int n, wchar_t *fa) {
		return call_dirsh(c, n, from_wchar(fa));
	}
}
