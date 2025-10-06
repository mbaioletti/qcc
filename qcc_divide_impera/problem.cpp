#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <sstream>
#include <set>
#include "problem.h"

using namespace std;

string Gate::str() {
    stringstream s;
    s << type << param << " q[" << to_string(qb1) << "]";
    if(arity==2) s << ", q[" << to_string(qb2)  << "]";
    return s.str();
}    

void Problem::load_architecture(string filename) {
	fstream f;
	f.open(filename, ios::in);
	if(f.fail()) {
		cerr << "Can't open " << filename << endl;
		exit(1);
	}
	f >> n_ph_qubits;
	connected=new bool[n_ph_qubits*n_ph_qubits];
	for(int i=0; i<n_ph_qubits; i++)
		for(int j=0; j<n_ph_qubits; j++)
			connected[i*n_ph_qubits+j]=false;
	int q1,q2;
	while(f >> q1 >> q2) {
		connections.push_back(make_pair(q1, q2));
		connected[q1*n_ph_qubits+q2]=connected[q2*n_ph_qubits+q1]=true;
	}
	f.close();
    cout << "Read " << filename  << " with " << n_ph_qubits << " qubits and " << connections.size() << " connections" << endl;
            
	//distances = compute_all_shortest_paths();
    compute_all_shortest_paths();
}

void Problem::load_qasm(string fname) {
    fstream f;
    f.open(fname, ios::in);
    if(f.fail()) {
        cerr << "Can't open " << fname << endl;
        exit(1);
    }
    //used_qubits=ArrayList<Int>()
    //var i=1
    regex gate_regex("\\s*(\\w+)\\s*(\\(.*\\))?\\s*q\\[(\\d+)\\](\\s*,\\s*q\\[(\\d+)\\])?;?\\s*");
    regex measure_regex("\\s*measure\\s*q\\[(\\d+)\\]\\s*->c\\[(\\d+)\\];?\\s*");
    num_gates=0;
    int gates_unary=0;
    int gates_binary=0;
    int last_qubit=-1;
    smatch match;
    int l=0;
    while(! f.eof()) {
        ++l;
        string line;
        getline(f, line);
        if(line=="") continue;
        stringstream input_stringstream(line);
        string a0;
        getline(input_stringstream, a0, ' ');
        if(a0=="OPENQASM" || a0=="include" || a0=="gate" || a0=="qreg" || a0=="creg" || line.substr(0,2)=="//")
            circuit_preamble += line +"\n";
        else {
            if(regex_match(line, match, measure_regex)) {
                string s1=match[0].str(), s2=match[1].str();
                measurements.push_back({stoi(s1), stoi(s2)});
            }
            else {
                if(regex_match(line, match, gate_regex)) {
                    string type=match[1].str();
                    string param=match[2].str();
                    int qb1=stoi(match[3].str());
                    string s3=match[5].str();
                    int qb2= s3=="" ? -1 : stoi(s3);
                    if(qb2>=0)
                        gates_binary++;
                    else
                        gates_unary++;
                    gates.push_back(Gate(num_gates,type,qb1,qb2,param));
                    num_gates++;
                    last_qubit=max(last_qubit, qb1);
                    last_qubit=max(last_qubit, qb2);
                }
                else {
                    // errore
                    cerr << "Errore in linea " << l << " :" << line << endl;
                    exit(1);
                }
            }
        }
    }
    n_log_qubits = last_qubit+1;
    compute_precedences();
    cout << "Read circuit " << fname << " with " << gates.size() << " gates (" << gates_unary << " unary and " << gates_binary << " binary), using " << n_log_qubits << " qubits, with depth " << depth << endl;
}

void Problem::compute_precedences() {
    //cout << "computing precedences" << endl;
    Gate** last=new Gate*[n_log_qubits];
    for(int i=0; i<n_log_qubits; i++) last[i]=nullptr;
    int* tm=new int[n_log_qubits];
    for(int i=0; i<n_log_qubits; i++) tm[i]=0;
    for(auto &g : gates) {
        Gate* p1=last[g.qb1];
        if(p1!=nullptr) {
            p1->successors.push_back(&g);
            g.predecessors.push_back(p1);
        }
        last[g.qb1]=&g;
        if(g.arity==2) {
            Gate* p2=last[g.qb2];
            if(p2!=nullptr) {
                p2->successors.push_back(&g);
                g.predecessors.push_back(p2);
            }
            last[g.qb2]=&g;
            int st = 1 + max(tm[g.qb1],tm[g.qb2]);
            tm[g.qb1] = tm[g.qb2]=st;
            g.level=st-1;
        }
        else {
            g.level=tm[g.qb1];
            tm[g.qb1]++;
        }
    }
    depth=tm[0];
    for(int i=1; i<n_log_qubits; i++)
        depth=max(depth, tm[i]);
    //cout << "size " << num_gates << " depth " << depth << endl;
    delete[] last;
    delete[] tm;
}



void Problem::breadth_first_search(int s) {
    vector<int> predecessor(n_ph_qubits);
    queue<int> q;
    for(int i=0; i<n_ph_qubits; i++)
        predecessor[i]=-1; // not yet visited
    predecessor[s]=s;
    set_distance(s,s,0);
    q.push(s);
    while(not q.empty()) {
        int u=q.front();
        q.pop();
        int d=get_distance(s,u);
        for(auto &conn : connections) 
            if(conn.first==u) {
                int w=conn.second;
                if(predecessor[w]<0) {
                    predecessor[w]=u;
                    set_distance(s, w, d+1);
                    q.push(w);
                }
            }
    }
    // reconstructs the paths from s to all the other nodes
    for(int i=0; i<n_ph_qubits; i++) 
        if(i==s) {
            int *path=new int[1];
            path[0]=s;
            set_path(s,i,path);
        }
        else {
            int d=get_distance(s,i);
            int *path=new int[d+1];
            int v=i;
            for(int k=d; k>=0; k--) {
                path[k]=v;
                v=predecessor[v];
            }
            set_path(s,i,path);
        }
}

void Problem::compute_all_shortest_paths() {
    distances=new int[n_ph_qubits*n_ph_qubits];
    paths=new int*[n_ph_qubits*n_ph_qubits];
    for(int s=0; s<n_ph_qubits; s++)
        breadth_first_search(s);
    /* to test distances and paths
    for(;;) {
        int u,v;
        cout << "insert two nodes ";
        cin >> u >> v;
        int d=get_distance(u,v);
        cout << "distance " << d << endl;
        cout << "shortest path ";
        int *p=get_path(u,v);
        for(int i=0; i<d+1; i++)
            cout << p[i] << " ";
        cout << endl;
    }*/
}

int* Problem::compute_distances() {
	int *m = new int[n_ph_qubits * n_ph_qubits];
	for(int i=0; i<n_ph_qubits * n_ph_qubits; i++)
		m[i]=0;
	for(auto conn : connections) {
        int x=conn.first, y=conn.second;
		m[x*n_ph_qubits+y] = 2;
		m[y*n_ph_qubits+x] = 2;
	}
	for(int i=0; i<n_ph_qubits; i++)
		m[i*n_ph_qubits+i]=1;
	bool changed=true;
	while (changed) {
		changed=false;
		for(int s=0; s<n_ph_qubits; s++) {
			for(int t=s+1; t<n_ph_qubits; t++) {
                for(auto conn : connections) {
                    int x=conn.first, y=conn.second;
					if (m[s*n_ph_qubits+x]>0 and m[y*n_ph_qubits+t]>0) {
						int d = m[s*n_ph_qubits+x] + m[y*n_ph_qubits+t];
						if (m[s*n_ph_qubits+t]==0 or d < m[s*n_ph_qubits+t]) {
							m[s*n_ph_qubits+t] = m[s*n_ph_qubits+x] + m[y*n_ph_qubits+t];
							m[t*n_ph_qubits+s] = m[t*n_ph_qubits+y] + m[x*n_ph_qubits+s];
							changed = true;
						}
					}
				}
			}
        }
    }
	return m;
}

/*
fun get_swap_id(x:Int, y:Int):Int {
    val res=swap_ids[x*num_swaps+y]
    if (res==-1) throw Throwable("Swap on ($x,$y) is not possible")
    return res
}
*/
	
