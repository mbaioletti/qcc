#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <string>
#include <sstream>

#include "solution.h"
using namespace std;
#include <filesystem>
namespace fs = filesystem;

struct {
    int max_num_swaps=INT_MAX;
    int prune=1;
    string obj_fun="swaps";
    double eps=0.1, alpha;
    string chooser="rwadd";
    string outfile="", resfile="";
    double timeout=10;
    Problem *problem;
    int check=0;
    int maxqubits=INT_MAX;
    int seed=-1;
    double delta=0.5;
    double beta=5, gamma=10;
    double beta1=1, beta2=1;
    double min_prob_reuse=0.1, max_prob_reuse=0.8;
    double p_try = 0.2;   
    double mut_rad = 1.0; 
    bool accurate=false;
    bool auto_adapt=false;
    string init_state="default";
    string tracefile="";       
    string instance_name="";   
    int num_sol_found=0;     
} compile_data;

mt19937 gen;
uniform_real_distribution<double> unif(0,1);

int a_lt_o=0, o_lt_a=0, a_eq_o=0;

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
        double best_dsum=values[0].dsum, best_dmin=values[0].dmin;
        for(int i=0; i<values.size(); i++) {
            if(values[i].dsum < best_dsum or (values[i].dsum == best_dsum and values[i].dmin < best_dmin)) {
                best_dsum=values[i].dsum;
                best_dmin=values[i].dmin;
                best.clear();
                best.push_back(i);
            } else if(values[i].dsum==best_dsum and values[i].dmin==best_dmin)
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
    

PActivity roulette_wheel_mult(Solution *s, vector<HActivity> &values) {
    int nv=values.size();
    int first_swap=-1;
    for(int i=0; i<nv; i++)
        if(values[i].cand->is_swap()) {
            first_swap = i;
            break;
    }
    // if there is some gates and some swaps, ignore all the swaps
    if(first_swap > 0) 
        nv = first_swap;
    /*cout << "Candidates " << endl;
    for(auto ha : values) {
        cout << ha.cand->str() << " " << ha.heur << endl;
    }
    int x; cin >> x;*/
    // normalization
    int max_dsum=values[0].dsum, 
        min_dsum=values[0].dsum, 
        max_dmin=values[0].dmin, 
        min_dmin=values[0].dmin;
    for(int i=1; i<nv; i++) {
        max_dsum=max(max_dsum, values[i].dsum);
        min_dsum=min(min_dsum, values[i].dsum);
        max_dmin=max(max_dmin, values[i].dmin);
        min_dmin=min(min_dmin, values[i].dmin);        
    }
    double sum_w=0;    
    vector<double> w(nv);
    double diff_dsum=max_dsum-min_dsum, diff_dmin=max_dmin-min_dmin;    
    for(int i=0; i<nv; i++) {
        //val tau=pheromone.value(v.cand,s)
        double dsum_n=diff_dsum==0.0 ? 0.5 : (values[i].dsum-min_dsum)/diff_dsum;
        double dmin_n=diff_dmin==0.0 ? 0.5 : (values[i].dmin-min_dmin)/diff_dmin;        
        w[i]=pow(1-dsum_n, compile_data.beta1) * pow(1-dmin_n, compile_data.beta2);
        //if(values[i].cand->is_swap())
        //    w[i] *= compile_data.delta;        
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

PActivity roulette_wheel_add(Solution *s, vector<HActivity> &values) {
    int nv=values.size();
    int first_swap=-1;
    for(int i=0; i<nv; i++)
        if(values[i].cand->is_swap()) {
            first_swap = i;
            break;
    }
    // if there is some gates and some swaps, ignore all the swaps
    if(first_swap > 0) 
        nv = first_swap;
    /*cout << "Candidates " << endl;
    for(auto ha : values) {
        cout << ha.cand->str() << " " << ha.heur << endl;
    }
    int x; cin >> x;*/
    // normalization
    int max_dsum=values[0].dsum, 
        min_dsum=values[0].dsum, 
        max_dmin=values[0].dmin, 
        min_dmin=values[0].dmin;
    for(int i=1; i<nv; i++) {
        max_dsum=max(max_dsum, values[i].dsum);
        min_dsum=min(min_dsum, values[i].dsum);
        max_dmin=max(max_dmin, values[i].dmin);
        min_dmin=min(min_dmin, values[i].dmin);        
    }
    double sum_w=0;    
    vector<double> w(nv);
    double diff_dsum=max_dsum-min_dsum, diff_dmin=max_dmin-min_dmin;    
    for(int i=0; i<nv; i++) {
        //val tau=pheromone.value(v.cand,s)
        double dsum_n=diff_dsum==0.0 ? 0.5 : (values[i].dsum-min_dsum)/diff_dsum;
        double dmin_n=diff_dmin==0.0 ? 0.5 : (values[i].dmin-min_dmin)/diff_dmin;        
        double eta=1-(compile_data.gamma * dsum_n + dmin_n)/(compile_data.gamma+1);
        //if(values[i].cand->is_swap())
        //    eta *= compile_data.swaps_discount;
        w[i]=pow(eta, compile_data.beta);
        //if(values[i].cand->is_swap())
        //    w[i] *= compile_data.delta;        
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

Solution *best = nullptr;
vector<Solution*> archive;

void insert_new_solution(Solution *s) {
    //for(int i=0; i<
}

Solution* find_solution(pair<int,int> *init_state, const chrono::steady_clock::time_point& deadline, double prob_reuse, double reuse)
{
    Solution *s;
    int tm = 0;
    if(best and unif(gen)<prob_reuse) {
        tm = int(best->makespan*reuse);
        s = Solution::partial_copy(best, tm);
    }
    else {
        s = new Solution(compile_data.problem, init_state);
        tm  = 0;
    }

    while (s->num_executed < compile_data.problem->num_gates &&
           s->num_swaps    <= compile_data.max_num_swaps)
    {
        if (chrono::steady_clock::now() >= deadline) {
            delete s;                 
            return nullptr;           
        }
        int ts = (compile_data.obj_fun == "swaps") ? INT_MAX : tm;
        auto values = s->evaluate_executable(ts, compile_data.prune, compile_data.accurate);

        if (!values.empty()) {
            PActivity act;
            if (compile_data.chooser == "greedy")  act = greedy_randomized(s, values);
            else if (compile_data.chooser == "rwmult")  act = roulette_wheel_mult(s, values);
            else act = roulette_wheel_add(s, values);
            s->execute(act);
        } else {
            ++tm;  
        }
    }
    if (s->num_executed == compile_data.problem->num_gates) {
        s->sort_act();
        return s;
    }
    return nullptr;
}

static bool file_exists(const string& path) {
   ifstream f(path);
    return f.good();
}

double log_uniform(double lo, double hi){
    //Logarithmic rescaling 
    uniform_int_distribution<> u(log10(lo), log10(hi));
    return pow(10.0, u(gen));
}

inline double clamp(double x, double lo, double hi){
     return max(lo, min(x, hi));
}

vector<int> parseNumbers(const string& input) {
    vector<int> numbers;
    stringstream ss(input);
    string item;

    while (getline(ss, item, ';')) {
        if (!item.empty()) {
            numbers.push_back(stoi(item));
        }
    }

    return numbers;
}

static Solution* optimize()
{
    Problem* p  = compile_data.problem;
    auto *init  = new pair<int,int>[p->n_ph_qubits];
    if(compile_data.init_state=="default") {
        for (int i = 0; i < p->n_ph_qubits; ++i) 
            init[i] = {i, i};
    } else {
        vector<int> is=parseNumbers(compile_data.init_state);
        int nq=is.size();
        for(int i=0; i<nq; i++)
            init[i] = { i, is[i] };
        for(int i=nq; i< p->n_ph_qubits; ++i) 
            init[i] = {i, i};
    }
    int duration = int(compile_data.timeout*1000);
    auto t_start  = chrono::steady_clock::now();
    auto deadline = t_start + chrono::milliseconds(duration);

    int    best_depth     = INT_MAX;
    int    best_swaps     = INT_MAX;
    double best_beta  = compile_data.beta;
    double best_gamma = compile_data.gamma;
    double best_delta = compile_data.delta;

    constexpr double BETA_MIN  = 0.1, BETA_MAX  = 5.0;
    constexpr double GAMMA_MIN = 0.1, GAMMA_MAX = 10.0;
    constexpr double DELTA_MIN = 0.001, DELTA_MAX = 0.01;
    uniform_real_distribution<double> mut(-compile_data.mut_rad, compile_data.mut_rad);

    /* INITIAL RANDOMISATION --------------------------------*/
    if (compile_data.auto_adapt and (compile_data.chooser == "rwadd" || compile_data.chooser == "rwmult")) {
        compile_data.beta  = log_uniform(BETA_MIN , BETA_MAX );
        compile_data.gamma = log_uniform(GAMMA_MIN, GAMMA_MAX);
        compile_data.delta = uniform_real_distribution<double>(DELTA_MIN, DELTA_MAX)(gen);
    }

    int num_sol_found = 0, sum_depth = 0, sum_nswaps = 0;
    double prob_reuse = 0, reuse = 0.5;
    while (chrono::steady_clock::now() < deadline) {

        double old_beta = compile_data.beta;
        double old_gamma = compile_data.gamma;
        double old_delta = compile_data.delta;
        bool mutated = false;

        auto perturb = [&](double x) { return x * (1.0 + mut(gen)); };
        if (compile_data.auto_adapt and unif(gen) < compile_data.p_try) {
            compile_data.beta  = clamp(perturb(compile_data.beta ), BETA_MIN , BETA_MAX );
            compile_data.gamma = clamp(perturb(compile_data.gamma), GAMMA_MIN, GAMMA_MAX);
            compile_data.delta = clamp(perturb(compile_data.delta), DELTA_MIN, DELTA_MAX);
            mutated = true;
        }
        compile_data.beta1 = compile_data.beta;
        compile_data.beta2 = compile_data.gamma;
        double remaining_time=chrono::duration_cast<chrono::milliseconds>(deadline - chrono::steady_clock::now()).count();
        double remaining_ratio =  int(10*remaining_time / duration) / 10.0;
        prob_reuse = max(compile_data.min_prob_reuse, min(compile_data.max_prob_reuse, 1-remaining_ratio));
        
        Solution *s = find_solution(init, deadline, prob_reuse, reuse);
                 
        if (chrono::steady_clock::now() >= deadline) {
            if (s) delete s;
            break;  // overall time budget exhausted
        }

        if (s) {
            bool better = false;
            sum_depth  += s->makespan;
            sum_nswaps += s->num_swaps;

            if (compile_data.obj_fun == "swaps" && s->num_swaps < best_swaps) {
                better     = true;
                best_swaps = s->num_swaps;
            } else if (compile_data.obj_fun == "depth" && s->makespan < best_depth) {
                better     = true;
                best_depth = s->makespan;
            }

            // time absolute for start
            auto elapsed = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - t_start).count();
            
            // ROW for solution found
            if (!compile_data.tracefile.empty()) {
                bool exists = file_exists(compile_data.tracefile);
                ofstream fout(compile_data.tracefile, ios::out | ios::app);
                if (!exists) {
                    fout << "instance,elapsed_time_ms,sol_idx\n";
                }
                fout << compile_data.instance_name << ","
                    << elapsed                  << ","
                    << num_sol_found            << "\n"; 
            }
            if (better) {
                if (best) delete best;
                best = s;
                best_beta  = compile_data.beta;
                best_gamma = compile_data.gamma;
                best_delta = compile_data.delta;

                cout << (compile_data.obj_fun == "swaps" ? "swaps " : "depth ")
                    << (compile_data.obj_fun == "swaps" ? best_swaps : best_depth)
                    << " after " << num_sol_found << " solutions and " << elapsed
                    << " ms (seed " << compile_data.seed << ")" << endl;
            } else {
                if (compile_data.auto_adapt and mutated) {
                    compile_data.beta  = old_beta;
                    compile_data.gamma = old_gamma;
                    compile_data.delta = old_delta;
                }
                delete s;
            }
            ++num_sol_found;
        } else {
            if (compile_data.auto_adapt and mutated) {
                compile_data.beta  = old_beta;
                compile_data.gamma = old_gamma;
                compile_data.delta = old_delta;
            }
        }
    }

    cout << num_sol_found << " solutions found" << endl;
    cout << "average depth " << double(sum_depth)/num_sol_found << ", average swaps " << double(sum_nswaps)/num_sol_found << endl;

    /* propagate best hyper‑parameters */
    compile_data.beta  = best_beta;
    compile_data.gamma = best_gamma;
    compile_data.delta = best_delta;

    if (!best || compile_data.check == 0) return best;

    string result = best->check();
    if (result != "ok") {
        cerr << "Compilation error: " << result << endl;
        best->save_qasm("compilation_error.qasm");
        return nullptr;
    }
    cout << "The solution found is correct" << endl;
    return best;
}


void read_options(int argc, char *argv[]);

string base_filename(string path) {
    return path.substr(path.find_last_of("/") + 1);
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: "
             << argv[0]
             << " <architecture_file> <qasm_file> [options]\n"
             << "  options:\n"
             << "    -objf <swaps|depth>\n"
             << "    -choose <greedy|rwadd|rwmult>\n"
             << "    -prune <int>\n"
             << "    -timeout <double>\n"
             << "    -seed <int>\n"
             << "    -res <result_csv>\n";
        return 1;
    }

    string arch_file = argv[1];
    string qasm_file = argv[2];

    read_options(argc - 3, argv + 3);

    if (compile_data.seed == -1) compile_data.seed = time(nullptr);
    cout << "Seed " << compile_data.seed << "\n";
    gen.seed(compile_data.seed);

    filesystem::path qp = qasm_file;
    compile_data.instance_name = qp.stem().string();

    Problem p;
    p.load_architecture(arch_file);
    p.load_qasm(qasm_file);
    compile_data.problem = &p;
    p.maxqubits = compile_data.maxqubits;

    Solution* s = optimize();

    if (compile_data.accurate) {
        cout << "accurate < original " << a_lt_o << ", accurate > original " << o_lt_a << ", accurate = original " << a_eq_o << endl;
    }

    if (s == nullptr) {
        cout << "No solution found" << endl;
    } else {
        cout << "Found a solution with " << s->activities.size() << " gates, depth " << s->makespan << " and " << s->num_swaps << " swaps " << endl;

        if (!compile_data.resfile.empty()) {
            bool exists = file_exists(compile_data.resfile);
            ofstream fout(compile_data.resfile, ios::out | ios::app);
            if (!exists) {
                fout << "instance_name,timeout,gate_out,depth_out,num_swaps,beta,gamma,delta,order,seed,num_sol_found\n";
            }
            fout
              << compile_data.instance_name        << ","
              << compile_data.timeout              << ","
              << s->activities.size()              << ","
              << s->makespan                       << ","
              << s->num_swaps                      << ","
              << compile_data.beta                 << ","
              << compile_data.gamma                << ","
              << compile_data.delta                << ","
              << compile_data.init_state           << ","
              << compile_data.seed                 << ","
              << compile_data.num_sol_found        << "\n";
            fout.close();
        }

        if (!compile_data.outfile.empty())
            s->save_qasm(compile_data.outfile);

        delete s;
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
        else if(opt=="-gamma") {
            compile_data.gamma=stod(val);
            i+=2;
        }                
        else if(opt=="-delta") {
            compile_data.delta=stod(val);
            i+=2;
        }                
        else if(opt=="-min_pr") {
            compile_data.min_prob_reuse=stod(val);
            i+=2;
        }                
        else if(opt=="-max_pr") {
            compile_data.max_prob_reuse=stod(val);
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
        else if(opt=="-accurate") {
            compile_data.accurate = stoi(val);
            i+=2;
        }      
        else if(opt=="-adapt") {
            compile_data.auto_adapt = stoi(val);
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
        else if(opt=="-init_state") {
            compile_data.init_state=val;
            i+=2;
        }
        else if(opt=="-trace") {
            compile_data.tracefile = val;
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
