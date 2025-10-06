#include <iostream>
#include <vector>
#include <unordered_set>
#include <fstream>
#include "problem.h"

//using namespace std;

enum ActivityState {
    NOT_EXECUTABLE=1,
    SUPPORTED=2,
    EXECUTABLE=3,
    EXECUTED=4
};

typedef struct Activity *PActivity;

struct Activity { 
    Gate *gate;
    ActivityState status = NOT_EXECUTABLE;
    int loc1=-1, loc2=-1;        
    int preconditions=0;
    int start_time=-1;
    int end_time=-1;
    int dur=0;
    PActivity just1= nullptr;
    PActivity just2= nullptr;
    Activity(Gate *g) {
       gate=g;
       if(gate->type != "swap" && gate->type != "start") 
            preconditions=gate->predecessors.size();
    }
    bool is_swap() { return gate->type == "swap"; }
    string str();
};

struct Solution {
    Problem *problem;
    pair<int,int> *init_state;
    int num_swaps=0;
    int num_executed=0;
    Gate *start_gate, *swap_gate;
    vector<PActivity> activities, possible_swaps;
    unordered_set<PActivity> not_executable, executable, supported;
    //vector<PActivity> zero_preconditions;
    int* assignment;
    int* position;
    PActivity* last_action;
    int makespan=0;
    int score=0;
    PActivity start;
    // methods
    Solution(Problem *p, pair<int,int> *is) {
        problem=p;
        init_state=is;
        allocate();
        init();
    }
    void init();
    void allocate();
    void deallocate();
    ~Solution() { deallocate(); }
    void create_initial_state();
    void concat(Solution&);
    void execute(PActivity);
    void save_qasm(string fname);
    void update_status();
    bool executable_at_time(PActivity, int);
    void evaluate_state_all_gates(int&, int&);
    void evaluate_state(int&, int&);
    vector<struct HActivity> evaluate_executable(int, int);
    vector<struct HActivity> prune(std::vector<HActivity>&, int, int);
    void swap_qubits(int l1,int l2);
    PActivity create_swap_activity(int lo1,int lo2);
    string check();
    void partial_copy(Solution *s,int last_time=INT_MAX);
};

struct HActivity {
    PActivity cand;
    int dsum, dmin;
};
