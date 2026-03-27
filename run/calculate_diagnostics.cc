#include "../code/include.hh"

#include <iostream>
#include <chrono>

using std::cout;
using std::endl;
using std::ofstream;

using namespace std::chrono;
int main(int argc, char* argv[]){

    int ls_bins = atoi(argv[2]);
    double eps_max = atof(argv[3]);
    
    linspace_and_gl* eps = new linspace_and_gl(0., eps_max, ls_bins, 5);

    int N_bins = ls_bins + 5;
    
    density* ics = new density(eps, IC_TCM, IC_NU_E, IC_NU_MU, IC_NUBAR_E, IC_NUBAR_MU, IC_MAX_DISTFUN);
    ics->set_T_Tcm(IC_TEMP, IC_TCM);

    collision_integral** integrators = new collision_integral*[N_bins];
    dep_vars* R = new dep_vars(N_bins);
    dep_vars* load = new dep_vars(N_bins);
    double results[2];
    
    for(int i = 0; i < N_bins; i++){
        if(atoi(argv[1]) == 0)
            integrators[i] = new nu_nu_collision(i, eps, true);
        else
            integrators[i] = new nu_e_collision(i, eps, true, IC_TCM);
    }
        
    for(int i = 0; i < N_bins; i++){
        integrators[i]->compute_R(IC_TCM, IC_TCM, results);
        R->set_value(i, max(results[0], results[1]));
        load->set_value(i, integrators[i]->estimate_load());
//        cout << i << ", " << R->get_value(i) << ", " << load->get_value(i) << endl;
    }

    ofstream file;
    file.open(argv[4]);

    R->print_csv(file);
    file << endl;
    eps->print_csv(file);
    file << endl;
    load->print_csv(file);
    file.close();

    for(int i = 0; i < N_bins; i++)
        delete integrators[i];
    delete[] integrators;
    
    delete R;
    
    delete eps;
    delete ics;


    return 0;
}