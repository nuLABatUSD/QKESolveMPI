#include "QKESolve.hh"
#include "density.hh"
#include "arrays.hh"
#include "constants.hh"

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

using std::max;
using std::ofstream;
using std::cout;
using std::endl;
using std::atof;
using std::string;


int main(int argc, char* argv[]){
    if(argc != 8){
        cout << "Appropriate usage ./program T_CM K_nue K_numu K_nubare K_nubarmu delta_m_squared file_header" << endl << "Abort" << endl;
        return 1;
    }

    linspace_and_gl* eps = new linspace_and_gl(0., 20., 201, 5);

    density* ic = new density(eps, atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), 0.9);
    
    QKE* sim = new QKE(eps, 0.8, atof(argv[6]), ic);
    
    string run_filename = string(argv[7]) + "_run.csv";
    string eps_filename = string(argv[7]) + "_eps.csv";
    
    sim->run(1000, 30, 5e20, run_filename, true);

    ofstream eps_file;
    eps_file.open(eps_filename);
    
    eps_file << "# ";
    for(int i = 1; i < 6; i++)
        eps_file << argv[i] << ", ";
    eps_file << argv[6] << endl;
    
    eps->print_csv(eps_file);
    
/*    for(int i = 0; i < eps->get_len()-1; i++)
        eps_file << eps->get_value(i) << ", ";
    eps_file << eps->get_value(eps->get_len()-1) << std::endl;
    
    for(int i = 0; i < eps->get_len()-1; i++)
        eps_file << eps->get_weight(i) << ", ";
    eps_file << eps->get_weight(eps->get_len()-1) << std::endl;*/
    
    eps_file.close();
    
    delete sim;
    delete ic;
    delete eps;

    return 0;
}

