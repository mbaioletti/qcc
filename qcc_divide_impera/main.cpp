#include "compiler.h"

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
                  << "    -res <result_csv> \n"
                  << "    -divide num_chunks \n"
                  << "    -init_state default|random|ils"
                  << endl;
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
          << compile_data.seed                    << ","              
          << compile_data.is_method               << ","
          << compile_data.obj_fun                 << endl;
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
                fout << "instance_name,timeout,num_chunks,gate_out,depth_out,num_swaps,seed,method_is,objfun" << endl;
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
              << compile_data.seed                    << ","
              << compile_data.is_method               << ","
              << compile_data.obj_fun                 << endl;
            fout.close();
        }
        if (not compile_data.outfile.empty())
            s->save_qasm(compile_data.outfile);
        delete s;
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
        else if(opt=="-seed_is") {
            compile_data.seed_init_state = stoi(val);
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
        else if(opt=="-init_state") {
            compile_data.is_method = val;
            i += 2;
        }
        else {
            cerr << "option " << opt << " not found" << endl;
            exit(1);
        }
    }
}
