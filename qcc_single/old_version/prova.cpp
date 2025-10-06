#include "solution.h"

int main(int argc, char *argv[]) {
    Problem p;

    p.load_architecture(argv[1]);    
    p.load_qasm(argv[2]);

    pair<int,int> *is=new pair<int,int>[p.n_ph_qubits];
    for(int i=0; i<p.n_ph_qubits; i++)
        is[i]=make_pair(i,i);
    Solution *s=new Solution(&p, is);
    s->save_qasm("prova.qasm");
}
