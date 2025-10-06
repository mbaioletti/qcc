#include <iostream>
#include <vector>
#include <climits>
#include <regex>
using namespace std;

struct Gate {
	int id;
	string type, param;
	int arity;
	int qb1, qb2;
    Gate(int _id, string _name, int _qb1, int _qb2, string _param="") {
        id=_id;
        type=_name;
        qb1=_qb1;
        qb2=_qb2;
        arity = type=="SWAP" || qb2>=0 ? 2 : 1;
        param=_param;
    }
	int level;
	vector<Gate*> predecessors, successors;
    string str();
};	

struct Measure {
    int qubit;
    int cbit;
};

struct Problem {
	vector<pair<int,int>> connections;
	bool* connected;
	int* distances;
    int** paths;
    vector<Gate> gates;
    int num_gates;
    int depth=0;
    int n_ph_qubits=0;
    int n_log_qubits=0;
    int maxqubits;
    string circuit_preamble="";
    vector<Measure> measurements;

    //lateinit var swap_ids:IntArray
    // methods    
	bool adjacent(int i, int j) {
		return connected[i*n_ph_qubits+j];
	}
    void load_architecture(string filename);
    void load_qasm(string filename);
    void compute_precedences();
    int* compute_distances();
    void compute_all_shortest_paths();
    void breadth_first_search(int s);
    void set_distance(int l1, int l2, int d) {
        distances[l1*n_ph_qubits+l2]=d;
    }
    int get_distance(int l1, int l2) {
        return distances[l1*n_ph_qubits+l2];
    }
    void set_path(int l1, int l2, int *path) {
        paths[l1*n_ph_qubits+l2]=path;
    }
    int* get_path(int l1, int l2) {
        return paths[l1*n_ph_qubits+l2];
    }    
    int binary_gate_cost(int l1, int l2) {
        return distances[l1*n_ph_qubits+l2]-1;
    }
};
