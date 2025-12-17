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

struct compiler_options {
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
    int seed_init_state=0;
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
};


Solution* optimize();

extern compiler_options compile_data;
extern int a_lt_o, o_lt_a, a_eq_o;
extern mt19937 gen;

