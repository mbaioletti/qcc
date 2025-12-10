#include "problem.h"
#include <fstream>

int main(int argc, char *argv[]) {    
    Problem p;
    string s = argv[1];
    //p.load_architecture(argv[1]);
    p.load_qasm(argv[1]);
    
    int i=s.rfind("/");
    
    fstream f("tutti.csv", ios::out | ios::app);
    f << s.substr(i+1) << " " << p.n_log_qubits << " " << p.gates.size() << " " << p.depth << endl;
    f.close();
}
