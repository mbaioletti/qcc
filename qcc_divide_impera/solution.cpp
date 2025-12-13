#include "solution.h"
#include <random>
#include <algorithm>
#include <unordered_map>
#include <cassert>
extern mt19937 gen;
extern uniform_real_distribution<double> unif;
extern int a_lt_o, o_lt_a, a_eq_o;

ostream& operator<<(ostream &s, const Activity &a) {
    if(a.gate->type=="swap") 
        s << "swap q[" << a.loc1 << "], q[" << a.loc2 << "]";
    else {
        s << a.gate->type << a.gate->param << " q[" << a.loc1 << "]";
        if(a.gate->arity==2) s << ", q[" << a.loc2 << "]";
    }
    return s;
}

string Activity::str() {
    stringstream s;
    if(gate->type=="swap") 
        s << "swap q[" << loc1 << "], q[" << loc2 << "]";
    else {
        s << gate->type << gate->param << " q[" << loc1<< "]";
        if(gate->arity==2) s << ", q[" << loc2 << "]";
    }
    return s.str();
}    
   
void Solution::init() {    
    remaining=0;
    for(auto &g : problem->gates) 
        if(g.level>=from_level and g.level<=to_level) {
            PActivity a=new Activity(&g);
            activities.push_back(a);
            from_gate_to_act[&g]=a;
            a->num_preconditions=0;
            ++remaining;
        }
    connect_activities();
    for(auto conn : problem->connections) {
        if(conn.first>=problem->maxqubits or conn.second>=problem->maxqubits)
            continue;
        bool found=false;
        for(auto sw : possible_swaps)
            if(sw->loc1==conn.second and sw->loc2==conn.first) {
                found=true;
                break;
            }
        if(found) continue;        
        Activity *a=new Activity(swap_gate);
        a->loc1=conn.first;
        a->loc2=conn.second;
        possible_swaps.push_back(a);
    }
    create_initial_state();
}

void Solution::connect_activities() {
    for(auto a : activities) {
        Gate* g1=a->gate;
        for(auto g2 : g1->successors) 
            if(g2->level<=to_level) {
                PActivity b=from_gate_to_act[g2];
                a->successors.push_back(b);
                ++b->num_preconditions;
            }
    }
    for(auto a : activities) {
        if(a->num_preconditions==0)
            supported.insert(a);
        else
            not_executable.insert(a);                
    }
}

void Solution::allocate() {
    assignment= new int [problem->n_ph_qubits];
    saved_assignment= new int [problem->n_ph_qubits];
    position  = new int [problem->n_ph_qubits];
    saved_position  = new int [problem->n_ph_qubits];
    last_action= new PActivity[problem->n_ph_qubits];
    swap_gate=new Gate(-1,"swap",-1,-1);
    start_gate=new Gate(-1,"start",-1,-1);
    start=new Activity(start_gate);
}

void Solution::deallocate() {
    delete[] assignment;
    delete[] saved_assignment;
    delete[] position;
    delete[] saved_position;
    delete[] last_action;
    for(int i=0; i<activities.size(); i++)
        delete activities[i];
    for(int i=0; i<possible_swaps.size(); i++)
        delete possible_swaps[i];
    delete swap_gate;
    delete start_gate;
    delete start;
}

void Solution::create_initial_state() {
    start->end_time=0;
    for(int i=0; i<problem->n_ph_qubits; i++) {
        assignment[i]=-1;
        last_action[i]=start;
    }
     for(int i=0; i<problem->n_ph_qubits; i++) {
        int q=init_state[i].first, l=init_state[i].second;
        assignment[l]=q;
        position[q]=l;
        last_action[l]=start;
    }
    int l=problem->n_log_qubits;
    for(int i=0; i<problem->n_ph_qubits; i++)
        if(assignment[i]==-1) {
            assignment[i]=l;
            position[l]=i;
            l++;
        }
    update_status();
}

void Solution::append(Solution *s2) {
    assignment=s2->assignment;
    position=s2->position;
    int new_makespan=makespan;
    for(auto a : s2->activities) {
        a->start_time += makespan;
        a->end_time   += makespan;
        if(a->end_time > new_makespan)
            new_makespan=a->end_time;
        activities.push_back(a);
    }
    makespan = new_makespan;
    num_swaps += s2->num_swaps;
}

void Solution::save_qasm(string fname) {
    sort(activities.begin(), activities.end(), 
        [](PActivity a,PActivity b) { return a->start_time < b->start_time; } );
    fstream f; 
    f.open(fname, ios::out);
    if(f.fail()) {
        cerr << "Can't open " << fname << endl;
        return;
    }
    f << problem->circuit_preamble << endl;    
    //f << "OPENQASM 2.0;" << endl;
    //f << "include \"qelib1.inc\";" << endl;
    //f << "qreg q[" << problem->n_ph_qubits << "];" << endl;
    //f << "qreg c[" << problem->n_log_qubits << "];" << endl;
    for(auto &a : activities) {
        f << a->str() << ";"  << " // " << a->start_time << endl;
    }
    for(auto m : problem->measurements) {
        int l=position[m.qubit];
        f << "measure q[" << l << "]->c[" << m.cbit << "];" << endl;
    }
    f.close(); 
}

/*
    fun save_final_state(fname:String) {
        val f=File(fname).printWriter()
        for(q in 0 until problem.n_log_qubits)
            f.println("${q} ${position[q]}")
        f.close()
    }
*/


void Solution::execute(PActivity a) {
    if(a->gate->type=="swap") {
        int lo1=a->loc1;
        int lo2=a->loc2;
        Gate* g=a->gate;
        g->qb1=assignment[lo1];
        g->qb2=assignment[lo2];
        a->just1=last_action[lo1];
        a->just2=last_action[lo2];
        a->start_time=max(a->just1->end_time, a->just2->end_time);
        a->dur=1;  //a->dur=problem.compute_duration(a->gate, a->loc1, a->loc2)
        last_action[lo1]=last_action[lo2]=a;
        swap_qubits(lo1,lo2);
        num_swaps++; }
    else {
        if (a->gate->arity==1) {
            a->loc1=position[a->gate->qb1];
            assert(a->loc1 >= 0);
            a->just1=last_action[a->loc1];
            a->start_time = a->just1->end_time;
            last_action[a->loc1]=a;
            a->dur=1; //a->dur=problem.compute_duration(a->gate, a->loc1)
        }
        else {
            a->loc1=position[a->gate->qb1];
            a->loc2=position[a->gate->qb2];
            assert(a->loc1 >= 0 and a->loc2 >= 0);
            a->just1=last_action[a->loc1];
            a->just2=last_action[a->loc2];
            a->start_time=max(a->just1->end_time, a->just2->end_time);
            a->dur=1;  //a->dur=problem.compute_duration(a->gate,a->loc1,a->loc2);
            last_action[a->loc1]=last_action[a->loc2]=a; 
        }
        ++num_executed;
        --remaining;
        for (auto b : a->successors) {
            --b->num_preconditions;
            if(b->num_preconditions==0) {
                not_executable.erase(b);
                supported.insert(b);
                b->status=SUPPORTED;
            }
        }
    }
    a->end_time = a->start_time + a->dur;
    a->status = EXECUTED;
    if(a->gate->type!="swap") {
        executable.erase(a);
        supported.erase(a);
    }
    //activities.push_back(a);
    makespan = max(makespan, a->end_time);
    update_status();
}

void Solution::update_status() {
    //for (auto a : unjustified) {
    //    if (a->status!=EXECUTED && a->preconditions==0) {
    for(PActivity a : supported) {
        if(a->gate->arity==1) {
            a->status=EXECUTABLE;
            executable.insert(a);
        }
        else {
            int l1=position[a->gate->qb1];
            int l2=position[a->gate->qb2];
            if(problem->adjacent(l1,l2)) {
                a->status=EXECUTABLE;
                executable.insert(a);
            }
            else {
                a->status=SUPPORTED;
                executable.erase(a);
            }
        }
    }
}

Solution* Solution::partial_copy(Solution *s, int last_time) {
    Solution *ns=new Solution(s->problem, s->init_state, s->from_level, s->to_level);
    
    unordered_map<Gate*, PActivity> from_gate_to_act;
    for(auto a: ns->activities)
        from_gate_to_act[a->gate]=a;
        
    for(auto a : s->activities)
        if(a->start_time<=last_time) {
            PActivity a1;
            if(a->is_swap()) 
                a1=ns->create_swap_activity(a->loc1, a->loc2);
            else 
                a1=from_gate_to_act[a->gate];
            ns->execute(a1);
        }
    return ns;
}

//fun candidates_at_time(unj:List<PActivity>, tm:Int) =
//    unj.filter { a -> a->status==EXECUTABLE && executable_at_time(a,tm) }

/* TODO
 * controllare che gli swap siano effettivamente eseguite al livello corretto 
 */
 
bool Solution::executable_at_time(PActivity a,int tm) {
    int st;
    if (a->gate->type=="swap") {
        auto j1=last_action[a->loc1];
        auto j2=last_action[a->loc2];
        st=max(j1->end_time, j2->end_time); 
    }
    else if(a->gate->arity==2) {
        int l1=position[a->gate->qb1];
        int l2=position[a->gate->qb2];
        auto j1=last_action[l1];
        auto j2=last_action[l2];
        st=max(j1->end_time, j2->end_time);
    }
    else {
        int l1=position[a->gate->qb1];
        auto j1=last_action[l1];
        st=j1->end_time; 
    }
    a->start_time=st;
    /*cout << a->str() << " is executable at time " << st << ", but we are at " << tm << endl;
    int x; cin >> x; */
    return st<=tm;
} 

/*
 * void Solution::evaluate_state_all_gates(int &dsum, int &dmin) {
    dsum=0;
    dmin=INT_MAX;
    for (auto a : unjustified) {
        if (a->gate->arity==2 // a->status!=EXECUTED
                ) {
            int l1=position[a->gate->qb1];
            int l2=position[a->gate->qb2];
            int d=problem->binary_gate_cost(l1,l2);
            dsum += d;
            dmin = min(d, dmin);
        }
    }
} */   

int Solution::estimate_swap_count(vector<PActivity> &bin_supported) {
    int min_swaps=INT_MAX;
    constexpr int NUM=10;
    for(int r=1; r<=NUM; r++) {
        shuffle(bin_supported.begin(), bin_supported.end(), gen);
        // save assignment and position arrays
        for(int i=0; i<problem->n_ph_qubits; i++) {
            saved_assignment[i] = assignment[i];
            saved_position[i] = position[i];
        }
        int nsw=0;
        for(auto a : bin_supported) {
            int l1=position[a->gate->qb1];
            int l2=position[a->gate->qb2];            
            int d=problem->get_distance(l1, l2);
            int *p=problem->get_path(l1, l2);
            // trova il punto k (a caso ?)
            int k;
            if(d%2==1)
                k=(d-1)/2;
            else
                k=unif(gen)<0.5 ? d/2 : d/2-1;
            for(int j=0; j<k; j++)
                swap_qubits(p[j], p[j+1]);
            for(int j=d; j>=k+2; j--)
                swap_qubits(p[j], p[j-1]);
            assert(problem->adjacent(position[a->gate->qb1], position[a->gate->qb2]));
            nsw += d-1;
        }
        min_swaps = min(min_swaps, nsw);
        // restore assignment and position arrays
        for(int i=0; i<problem->n_ph_qubits; i++) {
            assignment[i] = saved_assignment[i];
            position[i] = saved_position[i];
        }
    }
    return min_swaps;
}

PActivity Solution::find_swap(int l1, int l2) {
    for(auto sw : possible_swaps) {
        if((sw->loc1==l1 and sw->loc2==l2) or (sw->loc1==l2 and sw->loc2==l1))
            return sw;
    }
    cerr << "Error in find swap " << l1 << " " << l2 << endl;
    exit(1);
}

void Solution::evaluate_state(vector<PActivity> &bin_supported, int &dsum, int &dmin) {
    dsum=0;
    dmin=INT_MAX;
    
    if(bin_supported.empty()) return;
    for (auto a : bin_supported) {
        int l1=position[a->gate->qb1];
        int l2=position[a->gate->qb2];
        int d=problem->binary_gate_cost(l1,l2);
        /*int *p=problem->get_path(l1, l2);
        assert(p[0]==l1);
        assert(p[d+1]==l2);
        PActivity sw1=find_swap(l1,p[1]), sw2=find_swap(p[d],l2);
        helpful_swaps.insert(sw1);
        helpful_swaps.insert(sw2); */
        dsum += d;
        dmin = min(d,dmin);
    }
}

void Solution::evaluate_state_accurate(vector<PActivity> &bin_supported, int &dsum, int &dmin) {
    dsum=0;
    dmin=INT_MAX;
    
    if(bin_supported.empty()) return;
    dsum=estimate_swap_count(bin_supported);
    
    int dsum_o=0;
    for (auto a : bin_supported) {
        int l1=position[a->gate->qb1];
        int l2=position[a->gate->qb2];
        int d=problem->binary_gate_cost(l1,l2);
        dmin = min(d,dmin);
        dsum_o += d;
    }
    if(dsum != dsum_o) {
        if(dsum<dsum_o)
            a_lt_o ++;
        else
            o_lt_a ++;
    }
    else
        a_eq_o++;
    //int x; cin >> x;
}

bool Solution::is_useful_swap(PActivity sw) {
    int lo1=sw->loc1, lo2=sw->loc2;
    for (auto a : supported) if(a->gate->arity==2 and a->status!=EXECUTABLE) {
        int l1=position[a->gate->qb1];
        int l2=position[a->gate->qb2];
        int d=problem->get_distance(l1, l2);
        if(lo1==l1 and problem->get_distance(lo2, l2)<d) return true;
        if(lo2==l1 and problem->get_distance(lo1, l2)<d) return true;
        if(lo1==l2 and problem->get_distance(lo2, l1)<d) return true;
        if(lo2==l2 and problem->get_distance(lo1, l1)<d) return true;
    }
    return false;
}

vector<HActivity> Solution::evaluate_executable(int tm, int pruning, bool accurate) {
    vector<PActivity> bin_supported;
    unordered_set<int> useful_qubits;
    for (auto a : supported) if(a->gate->arity==2 and a->status!=EXECUTABLE) {
        bin_supported.push_back(a);
        int l1=position[a->gate->qb1];
        int l2=position[a->gate->qb2];
        useful_qubits.insert(l1);
        useful_qubits.insert(l2); 
    }
    //unordered_set<PActivity> helpful_swaps;    
    int init_dsum, init_dmin;
    if(accurate)        
        evaluate_state_accurate(bin_supported, init_dsum, init_dmin);
    else
        evaluate_state(bin_supported, init_dsum, init_dmin);
        
    vector<HActivity> values;
    
    for (auto e : supported) if(e->status==EXECUTABLE) {
        if (executable_at_time(e,tm)) {
            int dmin = 1;
            int dsum = e->gate->arity==2 ? init_dsum-1 : init_dsum;
            values.push_back({e,dsum,dmin});
       }
    }

    //unordered_set<PActivity> not_used;
    //for (auto sw : helpful_swaps) {    
    for (auto sw : possible_swaps) {
        int lo1=sw->loc1, lo2=sw->loc2;
        if(useful_qubits.find(lo1)==useful_qubits.end() and useful_qubits.find(lo2)==useful_qubits.end())
        //if(not is_useful_swap(sw))
            continue;
        if (executable_at_time(sw,tm)) {
            int dsum, dmin;
            /*  cout << "check swap " << lo1 << "," << lo2 << " executability at time " << tm <<endl;
            cout << act_to_str(*last_action[lo1]) << " " << last_action[lo1]->end_time << endl;
            cout << act_to_str(*last_action[lo2]) << " " << last_action[lo2]->end_time << endl; */
            swap_qubits(lo1,lo2);
            if(accurate)        
                evaluate_state_accurate(bin_supported, dsum, dmin);
            else
                evaluate_state(bin_supported, dsum, dmin);
            swap_qubits(lo1,lo2);

            values.push_back({sw,dsum,dmin});
        }
    }
    if(pruning==1) {
        return prune(values, init_dsum, init_dmin);
    }
//    else if(pruning==2)
//        return prune2(values,init_dsum,init_dmin)
    else
        return values;
}


vector<HActivity> Solution::prune(vector<HActivity> &values,int init_dsum,int init_dmin) {
    //var res=values.filter { v -> v.cand.gate.type!="swap" }
    //if(res.size!=0) return res as ArrayList<HPActivity>
    vector<HActivity> res;
    for(auto v : values)
        if(v.dsum<init_dsum || (v.dsum==init_dsum && v.dmin<init_dmin))
            res.push_back(v);
    if(res.size()!=0) return res;
    for(auto v : values)
        if(v.dmin<init_dmin)
            res.push_back(v);        
    return res;
}
 
 

void Solution::swap_qubits(int l1,int l2) {
    int qb1=assignment[l1];
    int qb2=assignment[l2];
    position[qb1]=l2;
    position[qb2]=l1;
    assignment[l2]=qb1;
    assignment[l1]=qb2;
}

PActivity Solution::create_swap_activity(int lo1,int lo2) {
    int new_id=activities.size();
    int qb1=assignment[lo1];
    int qb2=assignment[lo2];
    //assert(sw.arity==2 && qb1>=0 && qb2>=0)
    //PActivity a=new Activity(new Gate(new_id,"swap",qb1,qb2));
    PActivity a=new Activity(swap_gate);
    a->loc1=lo1;
    a->loc2=lo2;
    activities.push_back(a);
    return a;
}

/* codice Kotlin che cerca uno stato iniziale con la ILS
 * da tradurre in C++ dopo averlo adattato alle strutture dati

    fun iterated_local_search(start_time:Long) {
        init_state=InitState(problem.num_loc,{i -> Pair(i,i)})
        heur_init_state=select_initial_state(init_state)
        println("initial estimated number of swaps ${heur_init_state.first}")
        var elapsed_time = System.currentTimeMillis()-start_time
        n_ils_tries=0
        n_accepted=0
        //val max_time_ils=(max_time*time_ils).toInt()
        //while (elapsed_time<max_time_ils) {
        while(n_ils_tries<=max_ils_iter && elapsed_time<max_time) {
            heur_init_state=try_improve_initial_state(init_state,heur_init_state)
            elapsed_time = System.currentTimeMillis()-start_time
        }
        println("num. ILS iterations ${n_ils_tries}, num. success ${n_accepted} after ${elapsed_time/1000.0} sec.")
        println("improved estimated number of swaps ${heur_init_state.first}")    
    }
    
    fun select_initial_state(state:InitState):Pair<Int,Int> {
        val places=IntArray(problem.num_loc, { i -> i })
        places.shuffle(rand)
        //print("Initial state ")
        for(i in 0 until problem.num_loc) {
            state[i]=Pair(i,places[i])
        }
        //println()
        var res=Pair(0,0)
        if(use_ls_is) {
            val s0=Solution(problem,state)
            res=local_search_initial_state(state, s0)
        }
        //print("Improved state ")
        //for(i in 0 until problem.num_loc)
        //    print(" ${state[i].second}")
        //println()
        return res
    }


    fun perturb(curr:InitState):InitState {
        var nsw=(problem.num_loc*perturb_strength).toInt()
        val new_ord=InitState(problem.num_loc,{i -> curr[i] })
        // perturbation
        n_ils_tries++
        while(nsw>0) {
            val j1=(0 until problem.num_loc).random(rand)
            val j2=(0 until problem.num_loc).random(rand)
            if(j1!=j2) {
                val p1=new_ord[j1].second
                val p2=new_ord[j2].second
                new_ord[j1]=Pair(j1,p2)
                new_ord[j2]=Pair(j2,p1)
                nsw--
            }
        }
        return new_ord
    }

    fun try_improve_initial_state(curr:InitState, h:Pair<Int,Int>):Pair<Int,Int> {
        val new_ord=perturb(curr)
        // local search
        val s0=Solution(problem,new_ord)
        //val unj=s0.unjustified()
        val h1=local_search_initial_state(new_ord, s0)
        // if better, accept it
        if(h1.first<h.first || (h1.first==h.first && h1.second<h.second)) {
            for(i in 0 until problem.num_loc)
                curr[i]=new_ord[i]
            //print("After LS ")
            //for(i in 0 until problem.num_loc)
            //    print(" ${curr[i].second}")
            //println()
            n_accepted++
            perturb_strength=0.2
            return h1
        }
        //if(n_ils_tries%10==0) {
        //    perturb_strength=Math.min(0.8,perturb_strength*1.1)
        //}
        return h   
    }
   
    fun local_search_initial_state(res:InitState, s0:Solution):Pair<Int,Int> {
        var improved:Boolean
        var (current_dsum, current_dmin)=s0.evaluate_state_all_gates()
        do {
            improved=false
            var best_dsum=current_dsum
            var best_dmin=current_dmin
            var best_i=0
            var best_j=0
            //print("(${best_dsum},${best_dmin}) ")
            for(i in 0..problem.num_loc-1) {
                var q=res[i].second
                for(j in 0..problem.num_loc-1) {
                    var p=res[j].second
                    if(p!=q) {
                        s0.swap_qubits(q,p)
                        val (dsum,dmin)=s0.evaluate_state_all_gates()
                        s0.swap_qubits(q,p)
                        if(dsum<best_dsum || (dsum==best_dsum && dmin<best_dmin)) {
                            best_dsum=dsum
                            best_dmin=dmin
                            best_i=i
                            best_j=j
                        }
                    }
                }
            }
            if(best_dsum<current_dsum || (best_dsum==current_dsum && best_dmin<current_dmin)) {
                current_dsum=best_dsum
                current_dmin=best_dmin
                val pi=res[best_i].second
                val pj=res[best_j].second
                s0.swap_qubits(pi, pj)
                res[best_i]=Pair(best_i, pj)
                res[best_j]=Pair(best_j, pi)
                improved=true
            }
        } while(improved);
        //println()
        return Pair(current_dsum,current_dmin)
    }
*/


string Solution::check() {
    sort(activities.begin(), activities.end(), 
        [](PActivity a,PActivity b) { return a->start_time < b->start_time; } );
    unordered_map<Gate*, PActivity> from_gate_to_act;
    for(auto a : activities)
        if(! a->is_swap() && a->gate->type != "start")
            from_gate_to_act[a->gate]=a;
    int* available=new int[problem->n_ph_qubits];
    for(int i=0; i<problem->n_ph_qubits; i++) {
        int q=init_state[i].first, l=init_state[i].second;        
        assignment[l]=q;
        position[q]=l;        
        last_action[l]=start;
        available[l]=0;
    }
    for(auto a : activities) {
        string g;
        if(! a->is_swap()) {
            g=a->gate->str();
            for (auto b1 : a->gate->predecessors) {
                PActivity b=from_gate_to_act[b1];
                //if(!found)
                if(b==nullptr)
                    return b1->str()+" is not executed, but it should precede "+a->gate->str();
                if (a->start_time<b->end_time)
                    return g+" must be executed after "+b1->str();
                
            }
        }
        else {
            g=a->str();
        }
        int l1=a->loc1, l2=a->loc2;
        int t1=available[l1];
        int t2=a->gate->arity==2 ? available[l2] : 0;
        int t=max(t1, t2);             
        if (a->start_time<t)
            return g+" cannot be executed at time "+to_string(a->start_time);
        if (a->gate->arity==2 && !problem->adjacent(l1, l2))
            return g+" cannot be executed because of connectivity";
        if (a->gate->type!="swap" && position[a->gate->qb1]!=l1)
            return g+" cannot be executed because the log. qubit "+to_string(a->gate->qb1)+" is not in phys. qubit "+to_string(l1)+", it is in "+to_string(position[a->gate->qb1]);
        if (a->gate->arity==2 && a->gate->type!="swap" && position[a->gate->qb2]!=l2)
            return g+" cannot be executed because the log. qubit "+to_string(a->gate->qb2)+" is not in phys. qubit "+to_string(l2)+", it is in "+to_string(position[a->gate->qb2]);
        available[a->loc1]=a->end_time;
        if(a->gate->arity==2) available[a->loc2]=a->end_time;
        if (a->gate->type=="swap") {
            int q1=assignment[l1], q2=assignment[l2];
            assignment[l1]=q2;
            assignment[l2]=q1;
            position[q1]=l2;
            position[q2]=l1;
        }
    }
    return "ok";
}


