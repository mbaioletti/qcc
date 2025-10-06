#include "problem.h"

int main(int argc, char **argv) {
    string fname_a=argv[1], fname_c=argv[2];
    
    Problem p;
    
    p.load_qasm(fname_c);
    p.load_architecture(fname_a);
}
