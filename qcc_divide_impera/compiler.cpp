#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <deque>
#include <cassert>
#include "solution.h"
using namespace std;
#include <filesystem>
#include <algorithm>
namespace fs = filesystem;

struct {
    int max_num_swaps=INT_MAX;
    int prune=1;
    string obj_fun="swaps";
    double eps=0.1, alpha;
    string chooser="rwmult";
    string outfile="", resfile="";
    double timeout=10;
    Problem *problem;
    int check=1;
    int maxqubits=INT_MAX;
    int seed=-1;
    double delta=1;
    double beta=5, gamma=10;
    double min_prob_reuse=0, max_prob_reuse=0;
    double reuse=0.5;
    double p_try = 0.2;   
    double mut_rad = 1.0; 
    bool accurate=false;
    bool auto_adapt=false;
    bool bandit=true;
    int max_archive=1;
    bool any_layer=true;
    int num_chunks=1;
} compile_data;

struct ArmStats {
    int beta, gamma;
    int pulls = 0;
    double reward_sum = 0.0;
    double mean() const { return pulls ? reward_sum / pulls : 0.0; }
};

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
        //double eta=1-(compile_data.gamma * dsum + dmin)/(compile_data.gamma+1);
        //if(values[i].cand->is_swap())
        //    eta *= compile_data.swaps_discount;
        w[i]=pow(1-dsum_n, compile_data.beta) * pow(1-dmin_n, compile_data.gamma);
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

Solution *best_sol = nullptr;
deque<Solution*> archive;

// extract a random element from a deque
Solution* random_element(deque<Solution*> &a) {
    int n=a.size();
    uniform_int_distribution<int> distr(0,n-1);
    int k=distr(gen);
    auto it=a.begin();
    while(k>0) {
        --k;
        ++it;
    }
    assert(it!=a.end());
    return *it;
}

int num_copies = 0, num_improved_copies = 0;
bool copied;

Solution* find_solution(pair<int,int> *init_state, int from_level, int to_level,
    const chrono::steady_clock::time_point& deadline, double prob_reuse, double reuse) {

    Solution *s;
    int tm = 0;
    /*cout << "Initial state ";
    for(int i=0; i<compile_data.problem->n_ph_qubits; ++i) {
        cout << init_state[i].first << " " << init_state[i].second << ", ";
    }
    cout << endl;*/
    if(compile_data.max_archive and not archive.empty() and unif(gen)<prob_reuse) {
        Solution *s1=random_element(archive);
        tm = int(s1->makespan*reuse);
        s = Solution::partial_copy(s1, tm);
        ++num_copies;
        copied = true;
    }
    else {
        s = new Solution(compile_data.problem, init_state, from_level, to_level);
        tm  = 0;
        copied = false;
    }

    while (s->remaining > 0) {
        //cout << "remaining " << s->remaining << endl;
        if (chrono::steady_clock::now() >= deadline) {
            //cout << "timeout" << endl;
            delete s;                 
            return nullptr;           
        }
        int ts = (compile_data.obj_fun == "swaps" and compile_data.any_layer) ? INT_MAX : tm;
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
        //s->check_state();
    }
    return s;
}


double log_uniform(double lo, double hi){
    //Logarithmic rescaling 
    uniform_int_distribution<> u(log10(lo), log10(hi));
    return pow(10.0, u(gen));
}

inline double clamp(double x, double lo, double hi){
     return max(lo, min(x, hi));
}

/*------------------------------------------------------------------
* ================== MULTI-ARMED BANDIT (UCB) SECTION ============
* If the -bandit option is enabled, we construct a grid of values
* discrete for (beta, gamma). Each pair is an “arm” of the bandit.
* The UCB (Upper Confidence Bound) algorithm is then used to.
* select, at each iteration, which arm to try:
*
* UCB = mean_reward + sqrt( 2 * ln(total_pulls) / pulls_arm )
*
* Where:
* - mean_reward: average reward observed for the arm.
* - total_pulls : total number of total draws
* - pulls_arm : number of pulls of the current arm.
*
* The root term encourages exploration of unproven arms,
* while mean_reward encourages those that are already promising (exploitation).
*-----------------------------------------------------------------*/
vector<ArmStats> arms;
int total_pulls = 0;

Solution* optimize_chunk(int from_level, int to_level, pair<int,int> *init, int duration) {
    auto t_start  = chrono::steady_clock::now();
    auto deadline = t_start + chrono::milliseconds(duration);
    best_sol=nullptr;
    

    // beta and gamma current
    int best_depth = INT_MAX, best_swaps = INT_MAX;

    int best_beta  = static_cast<int>(compile_data.beta);
    int best_gamma = static_cast<int>(compile_data.gamma);





    int num_sol_found = 0;
    long sum_depth = 0, sum_nswaps = 0;
    double prob_reuse = 0.1, reuse = compile_data.reuse;

    while (chrono::steady_clock::now() < deadline)
    {   
        // Choose arms with ucb rule
        if (compile_data.bandit) {
            int chosen = -1;
            double best_ucb = -numeric_limits<double>::infinity();
            for (int k = 0; k < (int)arms.size(); ++k) {
                if (arms[k].pulls == 0) { chosen = k; break; }
                double ucb = arms[k].mean() + sqrt(2.0 * log(double(total_pulls)) / arms[k].pulls);
                if (ucb > best_ucb) { best_ucb = ucb; chosen = k; }
            }
            compile_data.beta  = arms[chosen].beta;
            compile_data.gamma = arms[chosen].gamma;
        }
        else if (compile_data.auto_adapt && unif(gen) < compile_data.p_try) {
            uniform_int_distribution<int> betad(1,10), gammad(1,10);
            compile_data.beta  = betad(gen);
            compile_data.gamma = gammad(gen);
        }
        
        double remaining_time = chrono::duration_cast<chrono::milliseconds>(deadline - chrono::steady_clock::now()).count();
        double remaining_ratio = int(100.0*remaining_time/duration)/100.0;
        prob_reuse = compile_data.min_prob_reuse + 
            (compile_data.max_prob_reuse -compile_data.min_prob_reuse)*(1-remaining_ratio);
            
        Solution* s = find_solution(init, from_level, to_level, deadline, prob_reuse, reuse);
        if (!s) break;   

        // negative reward for minimize 
        double reward = (compile_data.obj_fun == "swaps")
                          ? -double(s->num_swaps)
                          : -double(s->makespan);
        if (compile_data.bandit) {
            for (auto& a : arms)
                if (a.beta == (int)compile_data.beta && a.gamma == (int)compile_data.gamma) {
                    ++a.pulls;  
                    a.reward_sum += reward;  
                    break;
                }
            ++total_pulls;
        }

        sum_depth  += s->makespan;
        sum_nswaps += s->num_swaps;
        ++num_sol_found;

        bool better = (compile_data.obj_fun == "swaps")? (s->num_swaps < best_swaps): (s->makespan  < best_depth);

        if (better) {
            if(copied) ++num_improved_copies;
            //if (best_sol) delete best_sol;
            best_sol   = s;
            best_swaps = s->num_swaps;
            best_depth = s->makespan;
            best_beta  = (int)compile_data.beta;
            best_gamma = (int)compile_data.gamma;
            archive.push_back(s);
            if(archive.size()>compile_data.max_archive) {
                delete archive.front();
                archive.pop_front();
            }

            /*
            auto elapsed = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - t_start).count();
            cout << (compile_data.obj_fun=="swaps" ? "swaps " : "depth ")
                 << (compile_data.obj_fun=="swaps" ? best_swaps : best_depth)
                 << "  |  beta=" << best_beta << "  gamma=" << best_gamma << " pr_reuse " << prob_reuse 
                 << (copied ? " (C) " : " ")
                 << "  after " << num_sol_found
                 << " sols, " << elapsed << " ms (seed " << compile_data.seed << ")" << endl;*/
        } else {
            delete s;
        }
    }

    /*if (num_sol_found) {
        cout << num_sol_found << " solutions found, avg depth "
             << double(sum_depth)/num_sol_found
             << ", avg swaps " << double(sum_nswaps)/num_sol_found << '\n';
    }
    
    cout << "Number of copies " << num_copies  << " improved copies " << num_improved_copies << endl;
    */
    
    for(auto a : archive)
        if(a!=best_sol)
            delete a;
    archive.clear();
    
    if (best_sol) {
        compile_data.beta  = best_beta;
        compile_data.gamma = best_gamma;
    }
    return best_sol;    
}

Solution* optimize() {
 
    Problem* p  = compile_data.problem;
    auto* init  = new pair<int,int>[p->n_ph_qubits];
    for(int i = 0; i < p->n_ph_qubits; ++i) init[i] = {i, i};
    if (compile_data.bandit) {
        constexpr int BETA_MIN = 0, BETA_MAX = 5, BETA_BINS = BETA_MAX - BETA_MIN + 1;
        constexpr int GAMM_MIN = 0, GAMM_MAX = 5, GAMMA_BINS = GAMM_MAX - GAMM_MIN + 1;

        for (int i = 0; i < BETA_BINS;  ++i)
            for (int j = 0; j < GAMMA_BINS; ++j) {
                int b = BETA_MIN  + i * (BETA_MAX  - BETA_MIN ) / (BETA_BINS  - 1);
                int g = GAMM_MIN  + j * (GAMM_MAX  - GAMM_MIN ) / (GAMMA_BINS - 1);
                //if(b<g)
                    arms.push_back({b, g});
            }
        compile_data.auto_adapt = false;
    }    
    int total_duration = static_cast<int>(compile_data.timeout * 1000);
    int num_lev_per_chunk = p->depth / compile_data.num_chunks;
    int duration = total_duration / compile_data.num_chunks;
    int from_level = 0, to_level = num_lev_per_chunk;
    auto new_init=init;
    Solution *sol;
    for(int c=1; c<=compile_data.num_chunks; c++) {
        cout << "Chunk #" << c << " from L" << from_level << " to L" << to_level << " timeout " << duration << " ms" << endl;
        Solution *cs=optimize_chunk(from_level, to_level, new_init, duration);
        if(cs==nullptr)
            return cs;
        if(c==1)
            sol=cs;
        else 
            sol->append(cs);
        from_level = to_level + 1;
        if(c<compile_data.num_chunks-1)
            to_level += num_lev_per_chunk;
        else
            to_level = p->depth;
        // update init with the final state of cs
        new_init  = new pair<int,int>[p->n_ph_qubits];
        for(int i = 0; i < p->n_ph_qubits; ++i)
            new_init[i] = { i, cs->position[i] };            
        cout << "swaps " << cs->num_swaps << " depth " << cs->makespan << endl;
    }
    string result = sol->check();
    if (result != "ok") {
        cerr << "Compilation error: " << result << '\n';
        sol->save_qasm("compilation_error.qasm");
        return nullptr;
    }
    cout << "The solution found is correct\n";           
    //cout << "swaps " << sol->num_swaps << " depth "  << sol->makespan << endl;          
    return sol;
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
                  << "    -res <result_csv> \n";
        return 1;
    }

    string arch_file = argv[1];
    string qasm_file = argv[2];

    read_options(argc - 3, argv + 3);

    // Seed for mt19937
    if (compile_data.seed == -1) compile_data.seed = time(nullptr);
    cout << "Seed " << compile_data.seed << "\n";
    gen.seed(compile_data.seed);

    Problem p;
    p.load_architecture(arch_file);
    p.load_qasm(qasm_file);
    compile_data.problem = &p;
    p.maxqubits = compile_data.maxqubits;

    Solution* s = optimize();
    
    if(compile_data.accurate) {
        cout << "accurate < original " << a_lt_o << ", accurate > original " << o_lt_a << ", accurate = original " << a_eq_o << endl;
    }

    if(s==nullptr) {
        cout << "No solution found" << endl;
        filesystem::path p=qasm_file;
        ofstream fout(compile_data.resfile, ios::out | ios::app);
        fout 
          << p.stem().c_str()                     << ","
          << compile_data.timeout                 << ","
          << compile_data.num_chunks              << ","
          << INT_MAX                              << ","
          << INT_MAX                              << ","
          << INT_MAX                              << ","
          //<< compile_data.beta                    << ","
          //<< compile_data.gamma                   << ","
          //<< compile_data.delta                   << ","
          << compile_data.seed                   << endl;
        fout.close();        
    }
    else {
        cout << "Found a solution with " << s->activities.size() << " gates, depth " << s->makespan << " and " << s->num_swaps << " swaps " << endl;
        cout << "beta " << compile_data.beta << " gamma " << compile_data.gamma << endl;
        if (not compile_data.resfile.empty()) {
            fstream fin(compile_data.resfile, ios::in);
            bool exists=not fin.fail();
            fin.close();
            if(not exists) {
                ofstream fout(compile_data.resfile, ios::out);
                fout << "instance_name,timeout,num_chunks,gate_out,depth_out,num_swaps,seed" << endl;
            }
            filesystem::path p=qasm_file;
            ofstream fout(compile_data.resfile, ios::out | ios::app);
            fout 
              << p.stem().c_str()                     << ","
              << compile_data.timeout                 << ","
              << compile_data.num_chunks              << ","
              << s->activities.size()                 << ","
              << s->makespan                          << ","
              << s->num_swaps                         << ","
              //<< compile_data.beta                    << ","
              //<< compile_data.gamma                   << ","
              //<< compile_data.delta                   << ","
              << compile_data.seed                   << endl;
            fout.close();
        }
        if (not compile_data.outfile.empty())
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
        else if (opt == "-bandit") {        
            compile_data.bandit = stoi(val);
            i += 2;
        }                
        else if(opt=="-delta") {
            compile_data.delta=stod(val);
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
        else if(opt=="-min_pr") {
            compile_data.min_prob_reuse=stod(val);
            i+=2;
        }                
        else if(opt=="-max_pr") {
            compile_data.max_prob_reuse=stod(val);
            i+=2;
        }                                
        else if(opt=="-reuse") {
            compile_data.reuse=stod(val);
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
        else if(opt=="-arch") {
            compile_data.max_archive = stoi(val);
            i+=2;
        }
        else if(opt=="-divide") {
            compile_data.num_chunks = stoi(val);
            i += 2;
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
