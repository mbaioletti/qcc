#include "solution.h"
#include <algorithm>
#include <unordered_map>

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
    for(auto &g : problem->gates) {
        PActivity a=new Activity(&g);
        activities.push_back(a);
        if(a->preconditions==0)
            supported.insert(a);
        else
            not_executable.insert(a);        
    }

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

void Solution::allocate() {
    assignment= new int [problem->n_ph_qubits];
    position  = new int [problem->n_ph_qubits];
    last_action= new PActivity[problem->n_ph_qubits];
    swap_gate=new Gate(-1,"swap",-1,-1);
    start_gate=new Gate(-1,"start",-1,-1);
    start=new Activity(start_gate);
}

void Solution::deallocate() {
    delete[] assignment;
    delete[] position;
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

void Solution::concat(Solution &s2) {
    assignment=s2.assignment;
    position=s2.position;
    int new_makespan=makespan;
    for(auto a : s2.activities) {
        a->start_time += makespan;
        a->end_time += makespan;
        if(a->end_time > new_makespan)
            new_makespan=a->end_time;
        activities.push_back(a);
    }
    makespan = new_makespan;
    num_swaps += s2.num_swaps;
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
            a->just1=last_action[a->loc1];
            a->start_time = a->just1->end_time;
            last_action[a->loc1]=a;
            a->dur=1; //a->dur=problem.compute_duration(a->gate, a->loc1)
        }
        else {
            a->loc1=position[a->gate->qb1];
            a->loc2=position[a->gate->qb2];
            a->just1=last_action[a->loc1];
            a->just2=last_action[a->loc2];
            a->start_time=max(a->just1->end_time, a->just2->end_time);
            a->dur=1;  //a->dur=problem.compute_duration(a->gate,a->loc1,a->loc2);
            last_action[a->loc1]=last_action[a->loc2]=a; 
        }
        num_executed++;
        for (auto b : a->gate->successors) {
            auto b1=activities[b->id];
            --b1->preconditions;
            if(b1->preconditions==0) {
                not_executable.erase(b1);
                supported.insert(b1);
                b1->status=SUPPORTED;
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

void Solution::partial_copy(Solution *s, int last_time) {
    unordered_map<Gate*, PActivity> from_gate_to_act;
    for(auto a: s->activities)
        from_gate_to_act[a->gate]=a;
        
    for(auto a : s->activities)
        if(a->start_time<=last_time) {
            PActivity a1;
            if(a->is_swap()) 
                a1=create_swap_activity(a->loc1, a->loc2);
            else
                a1=from_gate_to_act[a->gate];
            execute(a1);
        }
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


void Solution::evaluate_state(int &dsum, int &dmin) {
    dsum=0;
    dmin=INT_MAX;

    for (auto a : supported) if(a->gate->arity==2) {
        int l1=position[a->gate->qb1];
        int l2=position[a->gate->qb2];
        int d=problem->binary_gate_cost(l1,l2);
        dsum += d;
        dmin = min(d,dmin);
    }
}


vector<HActivity> Solution::evaluate_executable(int tm, int pruning) {
    int init_dsum, init_dmin;
    evaluate_state(init_dsum, init_dmin);
  
    vector<HActivity> values;
    //val used_qubits=HashSet<Int>()
    for (auto e : supported) if(e->status==EXECUTABLE) {
        if (executable_at_time(e,tm)) {
            int dmin = 1;
            int dsum = e->gate->arity==2 ? init_dsum-1 : init_dsum;
            values.push_back({e,dsum,dmin});
            /*if(e->gate->arity==2) {
                int l1=position[e->gate->qb1];
                int l2=position[e->gate->qb2];
                used_qubits.add(l1);
                used_qubits.add(l2);
            }*/
       }
    }

    for (auto sw : possible_swaps) {
        //if(!(lo1 in used_qubits) && !(lo2 in used_qubits))
        //    continue
        int dsum, dmin;
        int lo1=sw->loc1;
        int lo2=sw->loc2;
        if (executable_at_time(sw,tm)) {
            /*  cout << "check swap " << lo1 << "," << lo2 << " executability at time " << tm <<endl;
            cout << act_to_str(*last_action[lo1]) << " " << last_action[lo1]->end_time << endl;
            cout << act_to_str(*last_action[lo2]) << " " << last_action[lo2]->end_time << endl; */
            swap_qubits(lo1,lo2);
            evaluate_state(dsum, dmin);
            swap_qubits(lo1,lo2);

            values.push_back({sw,dsum,dmin});
        }
    }
    //print("before pruning ${values.size}")
    //return values
    if(pruning==1) {
        return prune(values, init_dsum, init_dmin);
    }
//    else if(pruning==2)
//        return prune2(values,init_dsum,init_dmin)
    else
        return values;
}

/*
fun evaluate_swaps(unj:List<PActivity>, pruning:Int=1):ArrayList<HPActivity> {

    val (init_dsum,init_dmin)=evaluate_state(unj)

    var k1min=Int.MAX_VALUE 
    var k2min=Int.MAX_VALUE
    var k1max=0
    var k2max=0
    val values=ArrayList<HPActivity>()
    val used_qubits=HashSet<Int>()

    for ((lo1,lo2) in problem.all_swaps) {
        //if(!(lo1 in used_qubits) && !(lo2 in used_qubits))
        //    continue 
        val e=PActivity(a_swap)
        e.loc1=lo1;
        e.loc2=lo2;
        if (executable_at_time(e,Int.MAX_VALUE)) {
            swap_qubits(lo1,lo2)
            val (dsum,dmin)=evaluate_state(unj)
            swap_qubits(lo1,lo2);
            k1min=min(k1min,dsum);
            k1max=max(k1max,dsum);
            k2min=min(k2min,dmin);
            k2max=max(k2max,dmin);
            values.add(HPActivity(e,dsum,dmin,0.0))
        }
    }
    val diff1=k1max.toDouble()-k1min.toDouble()
    val diff2=k2max.toDouble()-k2min.toDouble()
    for (v in values) {
        val r1=if(diff1==0.0) 1.0 else (v.dsum-k1min)/diff1
        val r2=if(diff2==0.0) 1.0 else (v.dmin-k2min)/diff2
        v.heur=10*r1+r2
    }
    //return values
    if(pruning==1)
        return prune(values,init_dsum,init_dmin)
    else if(pruning==2)
        return prune2(values,init_dsum,init_dmin)
    else
        return values
}
*/

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
/*                bool found=false;
                PActivity b=nullptr;
                for(auto a1 : activities) 
                    if(a1->gate==b1) {
                        found=true;
                        b=a1;
                        break;
                    }
                //val b=activities.find { a1 -> a1.gate==b1 } 
*/
                PActivity b=from_gate_to_act[b1];
                //if(!found)
                if(b1==nullptr)
                    return b1->str()+" is not executed, bust it should precede "+a->gate->str();
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


/*
    val gg=activities.map { a -> a->gate }
    for(g in problem.gates) {
        if (!(g in gg))
            return "Required gate ${g} is missing"
    }
    activities.sortBy { a -> a->start_time }
    val available=IntArray(init_state.size)
    for ((q,l) in init_state) {
        //val q=x.first
        //val l=x.second
        assignment[l]=q
        position[q]=l
    }
    for (a in activities) {
        if (a->gate->type!="swap") {
            // check temporal constraints among PS and MIX
            for (b1 in a->gate->predecessors) {
                val b=activities.find { a1 -> a1.gate==b1 }
                if(b==null)
                    return "${b1} is not executed"
                if (a->start_time<b.end_time)
                    return "${a} must be executed after ${b}"
            }
        }
        val locs=if(a->gate->arity==2) arrayOf(a->loc1,a->loc2) else arrayOf(a->loc1)
        val t=locs.map{l -> available[l]}.maxOrNull() ?: 0
        if (a->start_time<t)
            return "${a} cannot be executed at time ${a->start_time}"
        if (a->gate->arity==2 && !problem.adjacent(a->loc1,a->loc2))
            return "${a} cannot be executed because of connectivity"
        if (a->gate->type!="swap" && position[a->gate->qb1]!=a->loc1)
            return "${a} cannot be executed because ${a->gate->qb1} contains a wrong qstate"
        if (a->gate->arity==2 && a->gate->type!="swap" && position[a->gate->qb2]!=a->loc2)
            return "${a} cannot be executed because ${a->gate->qb2} contains a wrong qstate"
        for(l in locs)
            available[l]=a->end_time
        if (a->gate->type=="swap") {
            val l1=a->loc1
            val l2=a->loc2
            val q1=assignment[l1]
            val q2=assignment[l2]
            assignment[l1]=q2
            assignment[l2]=q1
            position[q1]=l2
            position[q2]=l1
        }
    }
    return "ok"
}
*/
