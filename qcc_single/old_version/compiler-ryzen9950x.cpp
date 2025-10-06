#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include "solution.h"
using namespace std;

struct {
    int max_num_swaps=INT_MAX;
    int prune=1;
    string obj_fun="swaps";
    double eps=0.1, alpha, beta=3;
    string chooser="rw";
    string outfile="", resfile="";
    double timeout=10;
    Problem *problem;
    int check=0;
    int maxqubits=INT_MAX;
    int seed=-1;
    double swaps_discount=0.05;
} compile_data;

mt19937 gen;
uniform_real_distribution<double> unif(0,1);

PActivity greedy_randomized(Solution *s, vector<HActivity> &values) {
    vector<int> best;
    /*cout << "Candidati " << values.size() << " al tempo " << tm << endl;
    for(HActivity& ha : values) {
        PActivity a=ha.cand;
        cout << act_to_str(*a) << " " << ha.heur << endl;
    }
    int x; cin >> x;*/
    int j;
    if(unif(gen)<compile_data.eps) {
        uniform_int_distribution<int> distr(0, values.size()-1);
        j=distr(gen);
    }
    else {
        double min_h=12;
        for(int i=0; i<values.size(); i++) {
            if(values[i].heur<min_h) {
                min_h=values[i].heur;
                best.clear();
                best.push_back(i);
            } else if(values[i].heur<min_h)
                best.push_back(i);
        }
        uniform_int_distribution<int> distr(0, best.size()-1);
        int k=distr(gen);
        j=best[k];
    }
    PActivity e=values[j].cand;
    //cout << "Scelta " << e->str() << endl;
    //cin >> x;*/
    //cout << e->gate->type << " " << e->gate->qb1 << " " << e->gate->qb2 << endl;
    if (e->gate->type=="swap")
        return s->create_swap_activity(e->loc1, e->loc2);
    else
        return e;
}
    

PActivity roulette_wheel(Solution *s, vector<HActivity> &values) {
    int nv=values.size();
    /*cout << "Candidates " << endl;
    for(auto ha : values) {
        cout << ha.cand->str() << " " << ha.heur << endl;
    }
    int x; cin >> x;*/
    vector<double> w(nv);
    vector<int> best;
    double sum_w=0, max_w=-1, num_ex=0;
    for(int i=0; i<nv; i++) {
        //val tau=pheromone.value(v.cand,s)
        double eta=1.0-values[i].heur/11.0;
        if(values[i].cand->is_swap())
            eta *= compile_data.swaps_discount;
        w[i]=pow(eta, compile_data.beta);
        sum_w += w[i];
    }
    int j=0;
    double r=unif(gen)*sum_w;
    while(j<nv && r>w[j]) {
        r -= w[j];
        ++j;
    }
    PActivity e=values[j].cand;
    //cout << "Scelta " << e->str() << endl;
    //cout << e->gate->type << " " << e->gate->qb1 << " " << e->gate->qb2 << endl;
    if (e->gate->type=="swap")
        return s->create_swap_activity(e->loc1, e->loc2);
    else
        return e;
}

void print_solution(Solution *s) {
    cout << "----------------------------------" << endl;
    cout << "Preambolo " << s->problem->circuit_preamble << endl;
    cout << "----------------------------------" << endl;    
    for(auto a : s->activities) {
        cout << a->str() << " " << a->start_time << endl;
    }
    cout << "----------------------------------" << endl;
}

Solution* find_solution(pair<int,int> *init_state) {
    Solution* s=new Solution(compile_data.problem, init_state);
    int tm=0;
    int lev=0;
    while(s->num_executed<compile_data.problem->num_gates && s->num_swaps<=compile_data.max_num_swaps) {
        int tt=(compile_data.obj_fun=="swaps") ? INT_MAX : tm;
        auto values = s->evaluate_executable(tt, compile_data.prune);
        if (values.size()>0) {
            PActivity act;
            if(compile_data.chooser=="greedy") 
                act=greedy_randomized(s, values);
            else
                act=roulette_wheel(s, values);
            //cout << act->str() << " selected at time " << tm << endl;
            s->execute(act);
            lev = max(lev, act->gate->level);
        }
        else {
            //cout << "no action at time " << tm << endl;
            ++tm;
        }
    }
    if(s->num_executed==compile_data.problem->num_gates) {
            return s;
    }
    else {
        cout << "Incomplete solution " << endl;
        return nullptr;
    }
}

Solution* optimize() {
    Problem *p=compile_data.problem;
    pair<int,int> *is=new pair<int,int>[p->n_ph_qubits];
    for(int i=0; i<p->n_ph_qubits; i++)
        is[i]=make_pair(i,i);   
    auto t1=chrono::steady_clock::now();
    int best_depth=INT_MAX, best_num_swaps=INT_MAX;
    Solution *best=nullptr, *s;
    int num_sol_found=0;
    while(true) {
        auto t2=chrono::steady_clock::now();
        s=find_solution(is);
        auto diff_time=chrono::duration_cast<chrono::milliseconds>(t2-t1).count();
        if(diff_time > 1000*compile_data.timeout)
            break;
        if(s!=nullptr) {
            if(compile_data.obj_fun=="swaps" && s->num_swaps<best_num_swaps) {
                best_num_swaps=s->num_swaps;
                if(best != nullptr)
                    delete best;
                best=s;
                cout << "swaps " << best_num_swaps << " after " << num_sol_found << " solutions and " << diff_time << " millisecs " << endl;
            }
            else if(compile_data.obj_fun=="depth" && s->makespan<best_depth) {
                best_depth=s->makespan;
                if(best != nullptr)
                    delete best;                
                best=s;
                cout << "depth " << best_depth << " after " << num_sol_found << " solutions and " << diff_time << " millisecs " << endl;
            }
            else {
                delete s;
            }
            ++num_sol_found;
        }
    }
    cout << num_sol_found << " solutions found" << endl;
    if(best==nullptr || compile_data.check==0) return best;
    string result=best->check();
    if(result != "ok") {
        cout << "Compilation error: " << result << endl;
        best->save_qasm("compilation_error.qasm");
        best=nullptr;
    }
    else {
        cout << "The solution found is correct " << endl;
    }
    return best;
}

void read_options(int argc, char *argv[]);

string base_filename(string path) {
    return path.substr(path.find_last_of("/") + 1);
}

int main(int argc, char*argv[]) {
    if(argc<3) {
        cerr << "compile_qasm filename_arch filename_qasm [options]" << endl;
        exit(1);
    }
    string fname_a=argv[1], fname_c=argv[2];
    read_options(argc-3, argv+3);
    
    Problem p;
    p.load_qasm(fname_c);
    p.load_architecture(fname_a);
    compile_data.problem = &p;
    p.maxqubits=compile_data.maxqubits;
    p.swaps_discount=compile_data.swaps_discount;
    if(compile_data.seed == -1)
        compile_data.seed=time(0);
    cout << "Seed " << compile_data.seed << endl;
    gen.seed(compile_data.seed);
    Solution *s=optimize();
    if(s==nullptr) {
        cout << "No solution found" << endl;
    }
    else {
        cout << "Found a solution with " << s->activities.size() << " gates, depth " << s->makespan << " and " << s->num_swaps << " swaps " << endl;
        if(compile_data.outfile!="")
            s->save_qasm(compile_data.outfile);
    }
    if(compile_data.resfile!="") {
        fstream f;
        f.open(compile_data.resfile, ios::out|ios::app);
        f << base_filename(fname_c) << "," 
          << base_filename(fname_a) << "," 
          << compile_data.timeout << "," 
          << s->makespan << "," 
          << s->num_swaps 
          << endl;
    }
}

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
	for(int i=0; i<n; i++) {
		string s=from_wchar(c[i].name);
		if(c[i].num_qubits==2)
			cout << s << " " << c[i].qb1 << " " << c[i].qb2 << endl;
		else
			cout << s << " " << c[i].qb1  << endl;
	}
}

void read_options(int argc, char *argv[]) {
    int i=0;
    while(i<argc) {
        string opt=argv[i], val=argv[i+1];
        if(opt=="-eps") {
            compile_data.eps=stod(val);
            i+=2;
        }
        else if(opt=="-beta") {
            compile_data.beta=stod(val);
            i+=2;
        }        
        else if(opt=="-swaps_discount") {
            compile_data.swaps_discount=stod(val);
            i+=2;
        }                
        else if(opt=="-objf") {
            compile_data.obj_fun=val;
            i+=2;
        }
        else if(opt=="-choose") {
            compile_data.chooser=val;
            i+=2;
        }   
        else if(opt=="-out") {
            compile_data.outfile=val;
            i+=2;
        }   
        else if(opt=="-res") {
            compile_data.resfile=val;
            i+=2;
        }           
        else if(opt=="-check") {
            compile_data.check = stoi(val);
            i+=2;
        }                        
        else if(opt=="-prune") {
            compile_data.prune = stoi(val);
            i+=2;
        }   
        else if(opt=="-timeout") {
            compile_data.timeout = stoi(val);
            i+=2;
        }                
        else if(opt=="-maxqubits") {
            compile_data.maxqubits = stod(val);
            i+=2;
        }                   
        else if(opt=="-seed") {
            compile_data.seed = stoi(val);
            i+=2;
        }                                
        else {
            cerr << "option " << opt << " not found" << endl;
            exit(1);
        }
    }
}

extern "C" {
	void compile(Gate_from_qiskit c[], int n) {
		compile_orig(c, n);
	}
}

